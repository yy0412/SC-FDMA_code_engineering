function [nVar,info] = mynrNullRENoiseEstimate(rxGridSlots,nullIndSlots,excludedIndSlots)
%mynrNullRENoiseEstimate Estimate received-grid noise power on known null REs.
%   NVAR = mynrNullRENoiseEstimate(RXGRIDSLOTS,NULLINDSLOTS) pools the
%   samples selected by NULLINDSLOTS and returns mean(abs(Ynull).^2).
%   These indices must identify REs whose transmitted value is known to be
%   zero. No sample-mean removal or degrees-of-freedom correction is used.
%
%   An optional EXCLUDEDINDSLOTS cell array can contain pilot or data RE
%   indices. The function rejects any overlap with the selected null REs.

if nargin < 3
    excludedIndSlots = cell(size(rxGridSlots));
end
if ~iscell(rxGridSlots) || isempty(rxGridSlots) || ...
        ~iscell(nullIndSlots) || numel(nullIndSlots) ~= numel(rxGridSlots) || ...
        ~iscell(excludedIndSlots) || numel(excludedIndSlots) ~= numel(rxGridSlots)
    error('mynrNullRENoiseEstimate:SlotInputs', ...
        ['Received grids, null indices, and optional excluded indices must ' ...
        'be matching nonempty cell arrays.']);
end

nSlots = numel(rxGridSlots);
gridSize = size(rxGridSlots{1});
samplesPerSlot = numel(rxGridSlots{1});
nullIndices = zeros(0,1);
nullSamples = complex(zeros(0,1));
nullCounts = zeros(1,nSlots);
slotNoisePower = nan(1,nSlots);

for s = 1:nSlots
    grid = rxGridSlots{s};
    if ~isnumeric(grid) || ~isequal(size(grid),gridSize)
        error('mynrNullRENoiseEstimate:GridSizeMismatch', ...
            'All received slot grids must be numeric and have the same size.');
    end
    ind = double(nullIndSlots{s}(:));
    validateattributes(ind,{'numeric'}, ...
        {'real','finite','integer','>=',1,'<=',numel(grid)}, ...
        'mynrNullRENoiseEstimate','nullIndSlots');
    if numel(unique(ind)) ~= numel(ind)
        error('mynrNullRENoiseEstimate:DuplicateIndices', ...
            'Null RE indices must be unique within each slot.');
    end

    excluded = double(excludedIndSlots{s}(:));
    if ~isempty(excluded)
        validateattributes(excluded,{'numeric'}, ...
            {'real','finite','integer','>=',1,'<=',numel(grid)}, ...
            'mynrNullRENoiseEstimate','excludedIndSlots');
        if any(ismember(ind,excluded))
            error('mynrNullRENoiseEstimate:IndexOverlap', ...
                'Known null REs must not overlap excluded pilot or data REs.');
        end
    end

    yNull = double(grid(ind));
    if any(~isfinite(yNull))
        error('mynrNullRENoiseEstimate:NonfiniteSamples', ...
            'Received null-RE observations must be finite.');
    end
    nullCounts(s) = numel(ind);
    if ~isempty(yNull)
        slotNoisePower(s) = mean(abs(yNull).^2);
    end
    nullIndices = [nullIndices;ind+(s-1)*samplesPerSlot]; %#ok<AGROW>
    nullSamples = [nullSamples;yNull]; %#ok<AGROW>
end

if isempty(nullSamples)
    error('mynrNullRENoiseEstimate:NoNullREs', ...
        'At least one known zero-transmit RE is required for noise estimation.');
end

nVar = mean(abs(nullSamples).^2);
info = struct('Formula','mean(abs(Ynull).^2)', ...
    'Interpretation', ...
    'noise plus residual leakage/interference on known null REs', ...
    'NullIndices',nullIndices,'NullSamples',nullSamples, ...
    'NullRECount',numel(nullSamples), ...
    'NullRECountsPerSlot',nullCounts, ...
    'PerSlotNoiseVariance',slotNoisePower);
end
