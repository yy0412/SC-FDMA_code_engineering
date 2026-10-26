function [frequencyOffsetHz,coarseFrequencyOffsetHz, ...
        residualFrequencyOffsetHz,info] = mynrEstimateFrequencyFromLSPhase( ...
        preliminaryLS,dmrsTimeSamples,sampleRate,refinementSlotLag,searchRangeHz)
%mynrEstimateFrequencyFromLSPhase Joint CFO from inter-slot DMRS phases.
%   For each slot pair, sum raw same-subcarrier products before taking
%   their phase, as in the phase-difference derivation. The coarse
%   adjacent-slot and fine long-lag stages then minimize equal-coefficient
%   complex phase-residual sums. Each slot pair uses its own actual DMRS
%   time difference; no extra magnitude weights or effective time are used.

validateattributes(preliminaryLS,{'numeric'},{'2d'},mfilename,'preliminaryLS');
[~,numberOfSlots] = size(preliminaryLS);
validateattributes(dmrsTimeSamples,{'numeric'}, ...
    {'vector','real','finite','numel',numberOfSlots}, ...
    mfilename,'dmrsTimeSamples');
validateattributes(sampleRate,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename,'sampleRate');
times = double(dmrsTimeSamples(:))/double(sampleRate);
if any(diff(times) <= 0)
    error('mynrEstimateFrequencyFromLSPhase:TimeOrder', ...
        'DMRS observation times must be strictly increasing.');
end
if nargin < 4 || isempty(refinementSlotLag)
    refinementSlotLag = 8;
end
validateattributes(refinementSlotLag,{'numeric'}, ...
    {'scalar','real','finite','positive','integer'}, ...
    mfilename,'refinementSlotLag');

if numberOfSlots < 2
    frequencyOffsetHz = 0;
    coarseFrequencyOffsetHz = 0;
    residualFrequencyOffsetHz = 0;
    info = struct('CoarsePairCount',0,'FinePairCount',0, ...
        'Method','unweighted phase-difference cost');
    return;
end
refinementSlotLag = min(max(2,double(refinementSlotLag)),numberOfSlots-1);
[coarsePhasors,coarseDt] = localPairs(preliminaryLS,times,1);
if isempty(coarsePhasors)
    error('mynrEstimateFrequencyFromLSPhase:NoAdjacentPairs', ...
        'Adjacent slots must share at least one usable DMRS subcarrier.');
end

if nargin < 5 || isempty(searchRangeHz)
    limitHz = 0.5/max(coarseDt);
    searchRangeHz = [-limitHz limitHz];
end
validateattributes(searchRangeHz,{'numeric'}, ...
    {'real','finite','vector','numel',2},mfilename,'searchRangeHz');
searchRangeHz = double(searchRangeHz(:).');
if searchRangeHz(2) <= searchRangeHz(1)
    error('mynrEstimateFrequencyFromLSPhase:SearchRange', ...
        'searchRangeHz must have strictly increasing endpoints.');
end

coarseGrid = linspace(searchRangeHz(1),searchRangeHz(2),1601);
coarseCost = localCostGrid(coarsePhasors,coarseDt,coarseGrid);
[~,coarseMinimum] = min(coarseCost);
coarseFrequencyOffsetHz = localRefineMinimum( ...
    coarsePhasors,coarseDt,coarseGrid,coarseMinimum);

[finePhasors,fineDt] = localPairs( ...
    preliminaryLS,times,refinementSlotLag);
if isempty(finePhasors)
    residualFrequencyOffsetHz = 0;
    fineGrid = 0;
    fineCost = [];
else
    halfWidthHz = 0.5/max(fineDt);
    fineGrid = linspace(-halfWidthHz,halfWidthHz,1601);
    corrected = finePhasors.*exp( ...
        -1i*2*pi*coarseFrequencyOffsetHz*fineDt);
    fineCost = localCostGrid(corrected,fineDt,fineGrid);
    [~,fineMinimum] = min(fineCost);
    residualFrequencyOffsetHz = localRefineMinimum( ...
        corrected,fineDt,fineGrid,fineMinimum);
end
frequencyOffsetHz = coarseFrequencyOffsetHz+residualFrequencyOffsetHz;
info = struct('Method','unweighted phase-difference cost', ...
    'CoarsePairCount',numel(coarsePhasors), ...
    'FinePairCount',numel(finePhasors), ...
    'CoarseTimeDifferencesSeconds',coarseDt, ...
    'FineTimeDifferencesSeconds',fineDt, ...
    'CoarseGridHz',coarseGrid,'CoarseCost',coarseCost, ...
    'FineResidualGridHz',fineGrid,'FineCost',fineCost, ...
    'SearchRangeHz',searchRangeHz,'RefinementSlotLag',refinementSlotLag);
end

function [phasors,deltaT] = localPairs(ls,times,lag)
phasors = complex(zeros(0,1));
deltaT = zeros(0,1);
for s = 1:size(ls,2)-lag
    p = s+lag;
    usable = find(isfinite(ls(:,s)) & isfinite(ls(:,p)) ...
        & abs(ls(:,s)) > 0 & abs(ls(:,p)) > 0);
    if isempty(usable)
        continue;
    end
    combined = sum(ls(usable,p).*conj(ls(usable,s)));
    if abs(combined) == 0
        continue;
    end
    phasors(end+1,1) = combined/abs(combined); %#ok<AGROW>
    deltaT(end+1,1) = times(p)-times(s); %#ok<AGROW>
end
end

function cost = localCostGrid(phasors,deltaT,gridHz)
cost = zeros(size(gridHz));
for q = 1:numel(gridHz)
    prediction = exp(1i*2*pi*gridHz(q)*deltaT);
    cost(q) = sum(abs(phasors-prediction).^2);
end
end

function estimate = localRefineMinimum(phasors,deltaT,gridHz,index)
if isscalar(gridHz)
    estimate = gridHz;
    return;
end
first = max(1,index-1);
last = min(numel(gridHz),index+1);
objective = @(f) sum(abs(phasors-exp(1i*2*pi*f*deltaT)).^2);
estimate = fminbnd(objective,gridHz(first),gridHz(last));
end
