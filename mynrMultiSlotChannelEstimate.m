function [hEstGrids,nVar,info] = mynrMultiSlotChannelEstimate(carrier,rxGridSlots,refIndSlots,refSymSlots,varargin)
%mynrMultiSlotChannelEstimate Joint custom-K channel estimate over slots.
%
%   Slots are concatenated in the OFDM-symbol dimension and processed by
%   mynrChannelEstimate once. This keeps the estimator logic identical to
%   the single-slot custom-K estimator while allowing time averaging and
%   interpolation to use DM-RS observations from adjacent slots.

nSlots = numel(rxGridSlots);
if nSlots == 0
    error('mynrMultiSlotChannelEstimate:EmptyInput', ...
        'rxGridSlots must contain at least one slot.');
end

[K,L,R] = size(rxGridSlots{1});
for s = 1:nSlots
    slotGrid = rxGridSlots{s};
    if size(slotGrid,1) ~= K || size(slotGrid,2) ~= L || size(slotGrid,3) ~= R
        error('mynrMultiSlotChannelEstimate:GridSizeMismatch', ...
            'All received slot grids must have the same size.');
    end
end

rxGridAll = cat(2,rxGridSlots{:});
[refIndAll,refSymAll] = localConcatenateReferences(K,L,nSlots,refIndSlots,refSymSlots); % Concate the indices and symbols for multi-slot
[hEstAll,nVar,info] = mynrChannelEstimate(carrier,rxGridAll,refIndAll,refSymAll,varargin{:});

hEstGrids = cell(1,nSlots);
for s = 1:nSlots
    slotSymbols = (s-1)*L + (1:L);
    hEstGrids{s} = hEstAll(:,slotSymbols,:,:);
end

info.JointSlots = nSlots;
info.SymbolsEstimated = size(rxGridAll,2);
end

function [refIndAll,refSymAll] = localConcatenateReferences(K,L,nSlots,refIndSlots,refSymSlots)
P = 1;
for s = 1:nSlots
    ind = double(refIndSlots{s}(:));
    if ~isempty(ind)
        P = max(P,ceil(max(ind)/(K*L)));
    end
end

refIndAll = zeros(0,1);
refSymAll = complex(zeros(0,1,'like',refSymSlots{1}));

for s = 1:nSlots
    ind = double(refIndSlots{s}(:));
    sym = refSymSlots{s}(:);
    if numel(ind) ~= numel(sym)
        error('mynrMultiSlotChannelEstimate:UnequalSymIndCount', ...
            'REFIND and REFSYM for slot %d must contain the same number of elements.',s);
    end
    if isempty(ind)
        continue;
    end

    [k,l,p] = ind2sub([K L P],ind);
    l = l + (s-1)*L;
    indAll = sub2ind([K L*nSlots P],k,l,p);
    refIndAll = [refIndAll; indAll(:)]; %#ok<AGROW>
    refSymAll = [refSymAll; sym(:)]; %#ok<AGROW>
end
end
