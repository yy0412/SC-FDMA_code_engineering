function [cfoCRLBHz2,info] = mynrCRLB( ...
    Carrier,PUSCH,NSlot,nVar,dmrsPowerPerRE,knownAk)
%mynrCRLB CFO Cramer-Rao lower bound for the complex DMRS model.
%   With KNOWNak supplied, the only unknown parameter is CFO and absolute
%   DMRS times are used. Without KNOWNak, each DMRS subcarrier has an
%   unknown complex coefficient; time is then centered separately on every
%   subcarrier. Omitting KNOWNak preserves compatibility with Main.m.

if nargin < 5 || isempty(dmrsPowerPerRE)
    dmrsPowerPerRE = 1;
end
validateattributes(dmrsPowerPerRE,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'dmrsPowerPerRE');
validateattributes(nVar,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'nVar');

knownAkMode = nargin >= 6 && ~isempty(knownAk);
if knownAkMode
    if isscalar(knownAk)
        knownAk = repmat(knownAk,Carrier.NCustomSubcarriers,1);
    end
    validateattributes(knownAk,{'numeric'}, ...
        {'vector','numel',Carrier.NCustomSubcarriers},mfilename,'knownAk');
    knownAk = knownAk(:);
else
    knownAk = ones(Carrier.NCustomSubcarriers,1);
end

cfoObservationTimeSamplesBySlot = cell(NSlot,1);
cfoObservationSubcarriersBySlot = cell(NSlot,1);
cfoObservationPilotPowersBySlot = cell(NSlot,1);

slotStartSample = 0;

for nslot = 0:NSlot-1
    slotIdx = nslot+1;
    carrierCRLB = Carrier;
    carrierCRLB.NSlot = nslot;

    slotOFDMInfo = mynrOFDMInfo(carrierCRLB);

    dmrsIndicesCRLB = mynrPUSCHDMRSIndices(carrierCRLB,PUSCH);

    dmrsSymbolsCRLB = sqrt(dmrsPowerPerRE)* ...
        mynrPUSCHDMRS(carrierCRLB,PUSCH);

    [dmrsSubcarriers,dmrsSymbolNumbers] = ind2sub( [Carrier.NCustomSubcarriers Carrier.SymbolsPerSlot], double(dmrsIndicesCRLB));

    % Start position of every OFDM symbol relative to the slot
    symbolStartSamples = [0;cumsum(double(slotOFDMInfo.SymbolLengths(1:end-1)))];

    % Time at the center of the useful part of every DMRS symbol. Retain
    % integer samples so interval groups can be identified exactly.
    dmrsTimeSamplesThisSlot = slotStartSample + ...
        symbolStartSamples(dmrsSymbolNumbers) + ...
        double(slotOFDMInfo.CyclicPrefixLengths(dmrsSymbolNumbers)) + ...
        double(slotOFDMInfo.Nfft)/2;

    cfoObservationTimeSamplesBySlot{slotIdx} = dmrsTimeSamplesThisSlot(:);
    cfoObservationSubcarriersBySlot{slotIdx} = dmrsSubcarriers(:);
    cfoObservationPilotPowersBySlot{slotIdx} = abs(dmrsSymbolsCRLB(:)).^2;

    slotStartSample = slotStartSample + sum(double(slotOFDMInfo.SymbolLengths));
end

cfoObservationTimeSamples = vertcat(cfoObservationTimeSamplesBySlot{:});
cfoObservationTimes = cfoObservationTimeSamples/double(slotOFDMInfo.SampleRate);
cfoObservationSubcarriers = vertcat(cfoObservationSubcarriersBySlot{:});
cfoObservationPilotPowers = vertcat(cfoObservationPilotPowersBySlot{:});

observationChannelPowers = abs(knownAk(cfoObservationSubcarriers)).^2;

if knownAkMode
    % Ak (including its phase at t=0) is known, so CFO is the only unknown
    % parameter. The known-Ak CRB uses absolute times and must not center
    % them; changing the time origin would also change the known Ak phase.
    cfoEffectiveTimes = cfoObservationTimes;
    crlbModel = 'known-Ak';
else
    % Backward-compatible nuisance-Ak model used by Main.m. An unknown
    % initial phase exists on every DMRS subcarrier, so eliminate it by a
    % per-subcarrier Schur complement, equivalent to centering the times.
    cfoEffectiveTimes = zeros(size(cfoObservationTimes));

    for subcarrier = unique(cfoObservationSubcarriers).'
        thisSubcarrier = cfoObservationSubcarriers == subcarrier;
        thisTime = cfoObservationTimes(thisSubcarrier);
        thisPower = cfoObservationPilotPowers(thisSubcarrier) .* ...
            observationChannelPowers(thisSubcarrier);
        weightedMeanTime = sum(thisPower .* thisTime)/sum(thisPower);
        cfoEffectiveTimes(thisSubcarrier) = thisTime-weightedMeanTime;
    end
    crlbModel = 'unknown-Ak-nuisance';
end

% Exact, potentially nonuniform DMRS times and actual pilot/channel powers.
cfoDMRSTimeInformation = sum(cfoObservationPilotPowers .* ...
    observationChannelPowers .* cfoEffectiveTimes.^2);

% Fisher information for f in Hz from all complex DMRS observations.
cfoFisherInformation = 2*(2*pi)^2/nVar*cfoDMRSTimeInformation;

cfoCRLBHz2 = 1/cfoFisherInformation;

uniqueDMRSTimeSamples = unique(cfoObservationTimeSamples,'stable');
info = struct( ...
    'ObservationTimeSamples',cfoObservationTimeSamples, ...
    'ObservationTimesSeconds',cfoObservationTimes, ...
    'EffectiveObservationTimesSeconds',cfoEffectiveTimes, ...
    'AdjacentDMRSTimeIntervalsSamples',diff(uniqueDMRSTimeSamples), ...
    'AdjacentDMRSTimeIntervalsSeconds', ...
        diff(uniqueDMRSTimeSamples)/double(slotOFDMInfo.SampleRate), ...
    'SampleRate',double(slotOFDMInfo.SampleRate), ...
    'PilotPowers',cfoObservationPilotPowers, ...
    'ChannelPowers',observationChannelPowers, ...
    'KnownAk',knownAk, ...
    'Model',crlbModel, ...
    'DMRSPowerPerRE',double(dmrsPowerPerRE), ...
    'DMRSTimeInformation',cfoDMRSTimeInformation, ...
    'FisherInformationHz',cfoFisherInformation);
end
