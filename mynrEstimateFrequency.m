function [frequencyOffsetHz,info] = mynrEstimateFrequency( ...
    preliminaryLS,dmrsTimeSamples,sampleRate,knownAk,searchRangeHz)
%mynrEstimateFrequency Direct complex-domain LS estimate of CFO.
%   The observation model is
%
%       preliminaryLS(k,s) = knownAk(k)*exp(1i*2*pi*f*t(s)) + noise(k,s).
%
%   knownAk contains the known complex channel coefficient at the chosen
%   time origin. No conjugate products and no phase unwrapping are used.
%   The one-dimensional nonlinear LS objective is searched globally on a
%   deterministic grid and then refined locally with FMINBND.

validateattributes(sampleRate,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'sampleRate');

[nSubcarriers,nSlots] = size(preliminaryLS);
validateattributes(dmrsTimeSamples,{'numeric'}, ...
    {'real','vector','finite','numel',nSlots},mfilename,'dmrsTimeSamples');

if nargin < 4 || isempty(knownAk)
    knownAk = ones(nSubcarriers,1,'like',preliminaryLS);
end

if isscalar(knownAk)
    knownAk = repmat(knownAk,nSubcarriers,1);
end

validateattributes(knownAk,{'numeric'}, ...
    {'vector','numel',nSubcarriers},mfilename,'knownAk');
knownAk = knownAk(:);

dmrsTimes = double(dmrsTimeSamples(:).')/double(sampleRate);

if nargin < 5 || isempty(searchRangeHz)
    uniqueTimes = unique(dmrsTimes,'sorted');
    if numel(uniqueTimes) < 2
        error('mynrEstimateFrequency:InsufficientTimes', ...
            'At least two distinct DMRS times are required.');
    end
    maxUnambiguousHz = 0.5/min(diff(uniqueTimes));
    searchRangeHz = [-maxUnambiguousHz maxUnambiguousHz];
elseif isscalar(searchRangeHz)
    validateattributes(searchRangeHz,{'numeric'}, ...
        {'real','finite','positive'},mfilename,'searchRangeHz');
    searchRangeHz = [-double(searchRangeHz) double(searchRangeHz)];
else
    validateattributes(searchRangeHz,{'numeric'}, ...
        {'real','vector','numel',2,'finite','increasing'}, ...
        mfilename,'searchRangeHz');
    searchRangeHz = double(searchRangeHz(:).');
end

referenceGrid = repmat(knownAk,1,nSlots);
timeGrid = repmat(dmrsTimes,nSubcarriers,1);

valid = isfinite(real(preliminaryLS)) & isfinite(imag(preliminaryLS)) & ...
    isfinite(real(referenceGrid)) & isfinite(imag(referenceGrid)) & ...
    abs(referenceGrid) > eps;

observations = preliminaryLS(valid);
references = referenceGrid(valid);
observationTimes = timeGrid(valid);

if isempty(observations)
    error('mynrEstimateFrequency:NoValidObservations', ...
        'No valid LS observations have a known nonzero Ak value.');
end

% Direct complex LS cost. With equal DMRS power, the preliminary-LS noise
% variances are equal, so no additional weighting is required.
costFunction = @(frequencyHz) sum(abs(observations - references .* ...
    exp(1i*2*pi*frequencyHz*observationTimes)).^2);

% The cost is multimodal over a wide CFO interval. A deterministic grid
% locates the global basin; FMINBND then refines the same LS objective.
timeSpan = max(observationTimes)-min(observationTimes);
if timeSpan > 0
    targetGridSpacingHz = 1/(16*timeSpan);
else
    targetGridSpacingHz = diff(searchRangeHz)/256;
end
numGridPoints = max(257,ceil(diff(searchRangeHz)/targetGridSpacingHz)+1);
numGridPoints = min(numGridPoints,20001);
frequencyGridHz = linspace(searchRangeHz(1),searchRangeHz(2),numGridPoints);
gridCost = zeros(size(frequencyGridHz));

for gridIdx = 1:numGridPoints
    gridCost(gridIdx) = costFunction(frequencyGridHz(gridIdx));
end

[~,bestGridIdx] = min(gridCost);

if bestGridIdx == 1 || bestGridIdx == numGridPoints
    frequencyOffsetHz = frequencyGridHz(bestGridIdx);
else
    lowerHz = frequencyGridHz(bestGridIdx-1);
    upperHz = frequencyGridHz(bestGridIdx+1);
    options = optimset('Display','off','TolX',1e-9);
    frequencyOffsetHz = fminbnd(costFunction,lowerHz,upperHz,options);
end

info = struct( ...
    'AngularFrequencyRadPerSecond',2*pi*frequencyOffsetHz, ...
    'Cost',costFunction(frequencyOffsetHz), ...
    'SearchRangeHz',searchRangeHz, ...
    'GridSpacingHz',frequencyGridHz(2)-frequencyGridHz(1), ...
    'NumObservations',numel(observations), ...
    'KnownAk',knownAk, ...
    'ObservationTimesSeconds',observationTimes);
end
