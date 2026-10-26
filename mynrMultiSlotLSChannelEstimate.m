function [hEstGrids,nVar,info] = mynrMultiSlotLSChannelEstimate( ...
        carrier,rxGridSlots,refIndSlots,refSymSlots,varargin)
%mynrMultiSlotLSChannelEstimate Legacy raw-LS front end for Wiener interpolation.
%   This compatibility wrapper forms raw HLS=Y./X observations at the DMRS
%   REs, estimates received-grid noise power on explicitly known null REs,
%   and delegates all interpolation to mynrMultiSlotWienerInterpolation.
%
%   New processing chains can supply a different channel observation (for
%   example the fully averaged mynrMultiSlotChannelEstimate output) directly
%   to mynrMultiSlotWienerInterpolation instead of using this wrapper.

parser = inputParser;
parser.FunctionName = 'mynrMultiSlotLSChannelEstimate';
addParameter(parser,'NullIndSlots',{},@iscell);
addParameter(parser,'ResidualMean',[0;0],@isnumeric);
addParameter(parser,'ResidualCovariance',[],@isnumeric);
addParameter(parser,'SlotNumbers',[],@isnumeric);
addParameter(parser,'SlotStartSamples',[],@isnumeric);
addParameter(parser,'CyclicPrefixFraction',0.5,@isnumeric);
addParameter(parser,'UseCache',true,@(x) islogical(x) && isscalar(x));
parse(parser,varargin{:});
opts = parser.Results;

if ~iscell(rxGridSlots) || isempty(rxGridSlots) || ...
        ~iscell(refIndSlots) || ~iscell(refSymSlots) || ...
        numel(refIndSlots) ~= numel(rxGridSlots) || ...
        numel(refSymSlots) ~= numel(rxGridSlots)
    error('mynrMultiSlotLSChannelEstimate:SlotInputs', ...
        ['Received grids, reference indices, and symbols must be matching ' ...
        'nonempty cell arrays.']);
end

S = numel(rxGridSlots);
K = double(carrier.NCustomSubcarriers);
L = double(carrier.SymbolsPerSlot);
pilotLSGrids = cell(size(rxGridSlots));
for s = 1:S
    grid = rxGridSlots{s};
    if ~isnumeric(grid) || size(grid,1) ~= K || size(grid,2) ~= L || ...
            numel(grid) ~= K*L
        error('mynrMultiSlotLSChannelEstimate:GridSizeMismatch', ...
            'Every received grid must be a K-by-L SISO matrix.');
    end
    ind = double(refIndSlots{s}(:));
    sym = refSymSlots{s}(:);
    validateattributes(ind,{'numeric'}, ...
        {'real','finite','integer','>=',1,'<=',K*L}, ...
        'mynrMultiSlotLSChannelEstimate','refIndSlots');
    if numel(ind) ~= numel(sym) || numel(unique(ind)) ~= numel(ind)
        error('mynrMultiSlotLSChannelEstimate:PilotIndices', ...
            'Each slot needs unique SISO pilot indices and one symbol per index.');
    end
    if any(~isfinite(sym)) || any(abs(sym) == 0)
        error('mynrMultiSlotLSChannelEstimate:PilotSymbols', ...
            'Reference symbols must be finite and nonzero.');
    end
    y = grid(ind);
    if any(~isfinite(y))
        error('mynrMultiSlotLSChannelEstimate:ReceivedPilots', ...
            'Received pilot observations must be finite.');
    end
    pilotGrid = complex(zeros(K,L,'like',grid));
    pilotGrid(ind) = y.*conj(sym)./abs(sym).^2;
    pilotLSGrids{s} = pilotGrid;
end

[nVar,nullInfo] = mynrNullRENoiseEstimate( ...
    rxGridSlots,opts.NullIndSlots,refIndSlots);
[hEstGrids,info] = mynrMultiSlotWienerInterpolation( ...
    carrier,pilotLSGrids,refIndSlots,refSymSlots,nVar, ...
    'ResidualMean',opts.ResidualMean, ...
    'ResidualCovariance',opts.ResidualCovariance, ...
    'SlotNumbers',opts.SlotNumbers, ...
    'SlotStartSamples',opts.SlotStartSamples, ...
    'CyclicPrefixFraction',opts.CyclicPrefixFraction, ...
    'UseCache',opts.UseCache);

info.Method = ...
    'DMRS LS, null-RE noise power, frequency Wiener, time Wiener';
info.ObservationSource = 'raw LS estimates at DMRS REs';
info.NoiseEstimationFormula = nullInfo.Formula;
info.NoiseInterpretation = nullInfo.Interpretation;
info.NullIndices = nullInfo.NullIndices;
info.NullSamples = nullInfo.NullSamples;
info.NullRECount = nullInfo.NullRECount;
info.NullRECountsPerSlot = nullInfo.NullRECountsPerSlot;
info.PerSlotNoiseVariance = nullInfo.PerSlotNoiseVariance;
info.PilotLS = info.PilotObservations;
info.PilotOutputUsesRawLS = true;
end
