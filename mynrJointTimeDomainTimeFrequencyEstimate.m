function [timingOffsetSamples,frequencyOffsetHz,info] = ...
    mynrJointTimeDomainTimeFrequencyEstimate(rxWaveform,txReferenceWaveform, ...
    referenceSampleMask,nominalStartSample,sampleRate, ...
    timingSearchSamples,frequencySearchRangeHz)
%mynrJointTimeDomainTimeFrequencyEstimate Joint time-domain timing/CFO LS search.
%   The estimator minimizes, over integer timing offset D and CFO f,
%
%     J(D,f) = sum_{n in P} |y(n0+D+n)-x(n)*exp(j*2*pi*f*t)|^2,
%
%   where P contains time samples known to the receiver. The CFO phase uses
%   the absolute capture sample t=(n0+D+n-1)/sampleRate, matching
%   mynrApplyTimeFrequencyOffset. No FFT-domain LS division, conjugate
%   product, or phase unwrapping is used.

validateattributes(rxWaveform,{'numeric'},{'2d','nonempty'}, ...
    mfilename,'rxWaveform');
validateattributes(txReferenceWaveform,{'numeric'},{'2d','nonempty'}, ...
    mfilename,'txReferenceWaveform');

if size(rxWaveform,2) ~= size(txReferenceWaveform,2)
    error('mynrJointTimeFrequencyEstimate:AntennaMismatch', ...
        'RX and reference waveforms must have the same number of columns.');
end

referenceLength = size(txReferenceWaveform,1);

if islogical(referenceSampleMask)
    validateattributes(referenceSampleMask,{'logical'}, ...
        {'vector','numel',referenceLength},mfilename,'referenceSampleMask');
    referenceSampleIndices = find(referenceSampleMask(:));
else
    validateattributes(referenceSampleMask,{'numeric'}, ...
        {'real','vector','integer','positive','<=',referenceLength}, ...
        mfilename,'referenceSampleMask');
    referenceSampleIndices = unique(referenceSampleMask(:),'sorted');
end

if isempty(referenceSampleIndices)
    error('mynrJointTimeFrequencyEstimate:EmptyReference', ...
        'At least one known reference sample is required.');
end

validateattributes(nominalStartSample,{'numeric'}, ...
    {'real','scalar','integer','positive','finite'}, ...
    mfilename,'nominalStartSample');
validateattributes(sampleRate,{'numeric'}, ...
    {'real','scalar','positive','finite'},mfilename,'sampleRate');
validateattributes(timingSearchSamples,{'numeric'}, ...
    {'real','vector','integer','finite','nonempty'}, ...
    mfilename,'timingSearchSamples');
validateattributes(frequencySearchRangeHz,{'numeric'}, ...
    {'real','vector','numel',2,'finite','increasing'}, ...
    mfilename,'frequencySearchRangeHz');

timingSearchSamples = unique(double(timingSearchSamples(:).'),'sorted');
frequencySearchRangeHz = double(frequencySearchRangeHz(:).');
nominalStartSample = double(nominalStartSample);
sampleRate = double(sampleRate);

referenceSamples = txReferenceWaveform(referenceSampleIndices,:);
referenceEnergy = sum(abs(referenceSamples).^2,'all');
captureEnergy = sum(abs(rxWaveform).^2,'all');
if referenceEnergy <= eps
    error('mynrJointTimeFrequencyEstimate:ZeroReferenceEnergy', ...
        'The selected reference samples have zero energy.');
end

% Stage 1: a coarse grid resolves each main lobe with about four samples.
% Stage 2 below refines the best coarse basin continuously with FMINBND.
% Keeping these controls explicit avoids treating the coarse grid as the
% final frequency resolution.
coarseFrequencyOversampling = 4;
minimumCoarseGridPoints = 33;
referenceTimeSpan = double(referenceSampleIndices(end) - ...
    referenceSampleIndices(1))/sampleRate;
if referenceTimeSpan > 0
    targetGridSpacingHz = ...
        1/(coarseFrequencyOversampling*referenceTimeSpan);
else
    targetGridSpacingHz = ...
        diff(frequencySearchRangeHz)/(minimumCoarseGridPoints-1);
end

numFrequencyGridPoints = max(minimumCoarseGridPoints, ...
    ceil(diff(frequencySearchRangeHz)/targetGridSpacingHz)+1);
numFrequencyGridPoints = min(numFrequencyGridPoints,20001);
frequencyGridHz = linspace(frequencySearchRangeHz(1), ...
    frequencySearchRangeHz(2),numFrequencyGridPoints).';
frequencyGridSpacingHz = frequencyGridHz(2)-frequencyGridHz(1);

% CZT evaluates sum c(m)*exp(-j*2*pi*f*m/Fs) on the requested frequency
% grid. A separate phase factor restores the absolute capture sample index.
cztStart = exp(1i*2*pi*frequencyGridHz(1)/sampleRate);
cztRatio = exp(-1i*2*pi*frequencyGridSpacingHz/sampleRate);

numTimingCandidates = numel(timingSearchSamples);
costSurface = inf(numFrequencyGridPoints,numTimingCandidates);
refinedFrequencyHz = nan(1,numTimingCandidates);
refinedCost = inf(1,numTimingCandidates);

for timingIdx = 1:numTimingCandidates
    timingCandidate = timingSearchSamples(timingIdx);
    candidateStartSample = nominalStartSample+timingCandidate;
    candidateSampleIndices = candidateStartSample-1+referenceSampleIndices;

    if candidateSampleIndices(1) < 1 || ...
            candidateSampleIndices(end) > size(rxWaveform,1)
        continue;
    end

    receivedSamples = rxWaveform(candidateSampleIndices,:);
    selectedCrossProducts = sum(receivedSamples.*conj(referenceSamples),2);

    % Preserve the original sample positions between separated DMRS symbols.
    crossSequence = complex(zeros(referenceLength,1,'like',rxWaveform));
    crossSequence(referenceSampleIndices) = selectedCrossProducts;
    localFrequencyCorrelation = czt(crossSequence, ...
        numFrequencyGridPoints,cztRatio,cztStart);

    candidateStartZeroBased = candidateStartSample-1;
    absoluteTimePhase = exp(-1i*2*pi*frequencyGridHz * ...
        candidateStartZeroBased/sampleRate);
    coherentScore = real(localFrequencyCorrelation.*absoluteTimePhase);
    % This is the full-capture LS cost for a shifted, zero-padded known
    % training template. ||y||^2 and ||x||^2 are common to every D; using a
    % candidate-dependent window energy would bias timing toward low-energy
    % windows and would no longer represent the stated joint LS problem.
    costSurface(:,timingIdx) = captureEnergy+referenceEnergy-2*coherentScore;

    [~,bestGridIdx] = min(costSurface(:,timingIdx));

    % Refine the continuous CFO for this integer timing candidate using the
    % same time-domain LS objective, then compare all refined 2-D candidates.
    absoluteSampleIndices = double(candidateSampleIndices-1);
    frequencyCost = @(frequencyHz) captureEnergy+referenceEnergy - ...
        2*real(sum(selectedCrossProducts.*exp(-1i*2*pi*frequencyHz * ...
        absoluteSampleIndices/sampleRate)));

    if bestGridIdx == 1 || bestGridIdx == numFrequencyGridPoints
        refinedFrequencyHz(timingIdx) = frequencyGridHz(bestGridIdx);
    else
        lowerFrequencyHz = frequencyGridHz(bestGridIdx-1);
        upperFrequencyHz = frequencyGridHz(bestGridIdx+1);
        options = optimset('Display','off','TolX',1e-9);
        refinedFrequencyHz(timingIdx) = fminbnd(frequencyCost, ...
            lowerFrequencyHz,upperFrequencyHz,options);
    end
    refinedCost(timingIdx) = frequencyCost(refinedFrequencyHz(timingIdx));
end

[minimumCost,bestTimingIdx] = min(refinedCost);
if ~isfinite(minimumCost)
    error('mynrJointTimeFrequencyEstimate:SearchOutsideCapture', ...
        'No timing candidate keeps all reference samples inside the capture.');
end

timingOffsetSamples = timingSearchSamples(bestTimingIdx);
frequencyOffsetHz = refinedFrequencyHz(bestTimingIdx);

info = struct( ...
    'Cost',minimumCost, ...
    'TimingSearchSamples',timingSearchSamples, ...
    'FrequencySearchRangeHz',frequencySearchRangeHz, ...
    'FrequencyGridHz',frequencyGridHz, ...
    'FrequencyGridSpacingHz',frequencyGridSpacingHz, ...
    'CoarseFrequencyOversampling',coarseFrequencyOversampling, ...
    'CostSurface',costSurface, ...
    'RefinedFrequencyHzByTiming',refinedFrequencyHz, ...
    'RefinedCostByTiming',refinedCost, ...
    'ReferenceSampleIndices',referenceSampleIndices, ...
    'NumReferenceSamples',numel(referenceSampleIndices), ...
    'ReferenceEnergy',referenceEnergy);
end
