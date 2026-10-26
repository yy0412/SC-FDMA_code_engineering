% LS时频偏估计，LS全平均信道估计/LS信道估计/LMMSE信道估计(已知噪声)，空值噪声估计，维纳插值
clear;clc;
tic;

SNRIn = 3.4; % SNR range (dB)
NTxAnts = 1; % Number of transmit antennas
NRxAnts = 1; % Number of receive antennas
NSlot = 28;  % Number of slots used by one TB
ChannelEstimationMode = 'practical'; % 'perfect' or 'practical'
ChannelAveragingWindow = [0 0];

%% Carrier and PUSCH Configuration
Carrier = mynrCarrierConfig;  % Carrier resource grid configuration
Carrier.SubcarriersPerPRB = 4;
Carrier.NSizeGrid = 1;  % Bandwidth in number of resource blocks
Carrier.SubcarrierSpacing = 120; % kHz
Carrier.CyclicPrefix = 'Normal';
Carrier.NCellID = 0;  % Cell identity

% Time/frequency-offset configuration for the selectable joint LS methods.
% 2.5e3 Hz is 2.5 kHz, not 2500 kHz. The occupied bandwidth is 480 kHz.
ApplyTimeFrequencyOffsets = false;
EnableTimeFrequencyCompensation = false;
JointEstimatorDomain = 'time'; % 'timeLS' or 'PhaseDiff'; 
ActualTimingOffsetSamples = 2;
ActualFrequencyOffsetHz = 2.5e3;   
DMRSPowerPerRE = 2; % Two comb-2 DMRS REs preserve the four-RE OFDM-symbol energy
%KnownAk = ones(Carrier.NCustomSubcarriers,1); % Unity-channel CRB reference

PUSCH = nrPUSCHConfig;
PUSCH.PRBSet = 0:Carrier.NSizeGrid-1; % PUSCH PRB allocation
PUSCH.SymbolAllocation = [0,Carrier.SymbolsPerSlot]; % PUSCH symbol allocation in each slot
PUSCH.MappingType = 'A'; % PUSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
PUSCH.NID = Carrier.NCellID; % Scrambling identifiers
PUSCH.RNTI = 1;
PUSCH.TransformPrecoding = true; % Enable/disable transform precoding
PUSCH.NumLayers = 1;  % Number of PUSCH transmission layers
PUSCH.TransmissionScheme = 'nonCodebook'; % Transmission scheme ('nonCodebook','codebook')
PUSCH.NumAntennaPorts = 1; % Number of antenna ports for codebook based precoding
PUSCH.TPMI = 0;  % Precoding matrix indicator for codebook based precoding
PUSCH.Modulation = 'QPSK'; % 'pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'

% PUSCH DM-RS configuration
PUSCH.DMRS.DMRSTypeAPosition = 2; % Mapping type A only. First DM-RS symbol position (2,3)
PUSCH.DMRS.DMRSLength = 1;        % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
PUSCH.DMRS.DMRSAdditionalPosition = 0;  % Additional DM-RS symbol positions (max range 0...3)
PUSCH.DMRS.DMRSConfigurationType = 1;   % DM-RS configuration type (1,2)
PUSCH.DMRS.NumCDMGroupsWithoutData = 2; % Number of CDM groups without data
PUSCH.DMRS.NIDNSCID = 0;                % Scrambling identity (0...65535)
PUSCH.DMRS.NSCID = 0;                   % Scrambling initialization (0,1)
PUSCH.DMRS.NRSID = 0;                   % Scrambling ID for low-PAPR sequences (0...1007)
PUSCH.DMRS.GroupHopping = 0;            % Group hopping (0,1)
PUSCH.DMRS.SequenceHopping = 0;         % Sequence hopping (0,1)

% LDPC parameters
TargetCodeRate = 679/1024; % Code rate used to calculate transport block size
LDPCDecodingAlgorithm = 'Belief propagation'; % Available algorithms: 'Belief propagation', 'Layered belief propagation', 'Normalized min-sum', 'Offset min-sum'
MaximumLDPCIterationCount = 50;

% HARQ process and rate matching/TBS parameters
XOverhead = 0; % Set PUSCH rate matching overhead for TBS (Xoh)
rvSeq = 0; % Set up redundancy version (RV) sequence. HARQ disabled - single transmission with RV=0, no retransmissions

% Get information about the baseband waveform after OFDM modulation step
waveformInfo = mynrOFDMInfo(Carrier);
occupiedBandwidthHz = double(Carrier.NCustomSubcarriers) * double(Carrier.SubcarrierSpacing)*1e3;

if NTxAnts ~= 1 || NRxAnts ~= 1
    error('FrequencyOffset_Check currently implements joint estimators for the configured SISO link.');
end

% Engineering acquisition interval, not a strict ambiguity limit of the
% full waveform (which also contains within-symbol CFO information).
maxDMRSCFOSearchHz = 0.5*double(Carrier.SlotsPerSubframe)*1e3;%应该再除以DMRS间隔符号数
timingOffsetSearchSamples = 0:4;
frequencyOffsetSearchRangeHz = [-maxDMRSCFOSearchHz maxDMRSCFOSearchHz];

[~,puschInfo] = mynrPUSCHIndices(Carrier,PUSCH); % Get PUSCH RE indices and information containing the bit capacity, symbol capacity, and DMRS OFDM symbol locations.
GPerSlot = puschInfo.G; % Length of codeword per slot from the UL-SCH transport channel
GTotal = GPerSlot*NSlot; % Total length of codeword in a transport block
MRB = numel(puschInfo.PRBSet)*NSlot; % Total number of PRB-slot unit
trBlkSize = nrTBS(PUSCH.Modulation,PUSCH.NumLayers,MRB, puschInfo.NREPerPRB,TargetCodeRate,XOverhead);

fprintf('Subcarriers/PRB = %d, occupied bandwidth = %.1f kHz, NSlot = %d, Nfft = %d, TBS = %d\n',Carrier.SubcarriersPerPRB, ...
    occupiedBandwidthHz/1e3,NSlot,waveformInfo.Nfft,trBlkSize);

%% Processing Loop
BLER = zeros(2,numel(SNRIn));
BLER(1,:) = SNRIn;

for snrIdx = 1:numel(SNRIn)
    % Reset the random number generator so that each SNR point will experience the same noise realization
    rng('default');
    SNRdB = SNRIn(snrIdx);
    SNR = 10^(SNRdB/10);
    N0 = 1/sqrt(NRxAnts*double(waveformInfo.Nfft)*SNR);
    nVar = N0^2*double(waveformInfo.Nfft);

    % Apply Monte-Carlo simulation for computing BLER with randomness
    MC = 10;
    blkerr = zeros(MC,1);
    
    % Log
    estimatedTimingOffsets = zeros(MC,1);
    estimatedFrequencyOffsetsHz = zeros(MC,1);
    TBSlog = zeros(trBlkSize,MC);
    codedTBSlog = zeros(GTotal,MC);
    
    %
    for mc = 1:MC
       % Create UL-SCH encoder System object to perform transport channel encoding
        encodeULSCH = nrULSCH;
        encodeULSCH.TargetCodeRate = TargetCodeRate;

        % Create UL-SCH decoder System object to perform transport channel decoding
        % Use layered belief propagation for LDPC decoding, with half the number of
        % iterations as compared to the default for belief propagation decoding
        decodeULSCH = nrULSCHDecoder;
        decodeULSCH.TargetCodeRate = TargetCodeRate;
        decodeULSCH.LDPCDecodingAlgorithm = LDPCDecodingAlgorithm;
        decodeULSCH.MaximumLDPCIterationCount = MaximumLDPCIterationCount;

        % Generate the transport block and coded block
        trBlk = randi([0 1],trBlkSize,1);
        setTransportBlock(encodeULSCH,trBlk);
        codedTrBlock = encodeULSCH(PUSCH.Modulation,PUSCH.NumLayers,GTotal,rvSeq);

        % Generate all transmit slots and save their DMRS information.
        txWaveformSlots = cell(1,NSlot);
        slotLengths = zeros(1,NSlot);
        puschIndicesSlots = cell(1,NSlot);
        puschInfoSlots = cell(1,NSlot);
        dmrsIndicesSlots = cell(1,NSlot);
        dmrsSymbolsSlots = cell(1,NSlot);
        
        TBSlog(:,mc) = trBlk;
        codedTBSlog(:,mc) = codedTrBlock;

        for nslot = 0:NSlot-1

            % Local parameters
            pusch = PUSCH;
            carrier = Carrier;
            carrier.NSlot = nslot;
            slotIdx = nslot+1;

            % Allocate the coded TB to each slot     
            slotRange = nslot*GPerSlot+(1:GPerSlot);
            codedTrBlockSlot = codedTrBlock(slotRange);
            
            % Generate the DMRS and its indices
            dmrsIndicesSlot = mynrPUSCHDMRSIndices(carrier,pusch);
            % Power boost: two occupied DMRS REs each have power 2, so the
            % DMRS OFDM symbol has the same total energy (4) as a data symbol.
            dmrsSymbolsSlot = sqrt(DMRSPowerPerRE)*mynrPUSCHDMRS(carrier,pusch);

            % Generate PUSCH symbol and allocate all RE (data+DMRS)
            puschGrid = mynrResourceGrid(carrier,NTxAnts);
            [puschIndicesSlot,puschInfoSlot] = mynrPUSCHIndices(carrier,pusch);
            puschGrid(puschIndicesSlot) = mynrPUSCH(carrier,pusch,codedTrBlockSlot);
            puschGrid(dmrsIndicesSlot) = puschGrid(dmrsIndicesSlot)+dmrsSymbolsSlot;

            txWaveformSlots{slotIdx} = mynrOFDMModulate(carrier,puschGrid);
            slotLengths(slotIdx) = size(txWaveformSlots{slotIdx},1);
            puschIndicesSlots{slotIdx} = puschIndicesSlot;
            puschInfoSlots{slotIdx} = puschInfoSlot;
            dmrsIndicesSlots{slotIdx} = dmrsIndicesSlot;
            dmrsSymbolsSlots{slotIdx} = dmrsSymbolsSlot;
        end

        %% Apply timing offset, frequency offset and AWGN
        % 添加时频偏
        slotStartSamples = [0 cumsum(slotLengths(1:end-1))];
        txWaveformTB = vertcat(txWaveformSlots{:});
        tbWaveformLength = size(txWaveformTB,1);

        % Preserve the existing time-estimator mask. The frequency branch
        % below constructs its reference from the known pilots and excludes
        % TX-window overlap with unknown neighboring OFDM symbols.
        dmrsTimeSampleMask = false(tbWaveformLength,1);
        dmrsUsefulSampleIndices = zeros(waveformInfo.Nfft,NSlot);
        % Known pilots only, in ordinary (unshifted) FFT-bin order. Do not
        % copy the windowed TX useful samples: their final Windowing samples
        % overlap the following unknown data symbol.
        dmrsFrequencySymbols = complex(zeros(waveformInfo.Nfft,NSlot));
        dmrsCyclicPrefixLengths = zeros(1,NSlot);
        dmrsWindowingSamples = zeros(1,NSlot);
        for nslot = 0:NSlot-1
            carrier = Carrier;
            carrier.NSlot = nslot;
            slotIdx = nslot+1;
            slotOFDMInfo = mynrOFDMInfo(carrier);

            [dmrsSubcarriers,dmrsSymbolNumbers] = ind2sub([Carrier.NCustomSubcarriers Carrier.SymbolsPerSlot], dmrsIndicesSlots{slotIdx});
            dmrsSymbolNumbers = unique(dmrsSymbolNumbers);
            if numel(dmrsSymbolNumbers) ~= 1
                error(['The joint estimator currently expects one DMRS ' 'OFDM symbol per slot.']);
            end
            dmrsSymbolNumber = dmrsSymbolNumbers(1);

            [~,dataSymbolNumbers] = ind2sub([Carrier.NCustomSubcarriers Carrier.SymbolsPerSlot],puschIndicesSlots{slotIdx});
            if any(dataSymbolNumbers == dmrsSymbolNumber)
                error(['The selected DMRS OFDM symbol contains unknown ' 'data RE and cannot be used as a fully known ' 'time-domain reference.']);
            end

            symbolStartInSlot = sum(double(slotOFDMInfo.SymbolLengths(1:dmrsSymbolNumber-1)));
            dmrsTimeRange = slotStartSamples(slotIdx)+symbolStartInSlot + (1:double(slotOFDMInfo.SymbolLengths(dmrsSymbolNumber)));
            dmrsTimeSampleMask(dmrsTimeRange) = true;
            usefulStartSample = dmrsTimeRange(1) + double(slotOFDMInfo.CyclicPrefixLengths(dmrsSymbolNumber));
            dmrsUsefulSampleIndices(:,slotIdx) = usefulStartSample + (0:double(slotOFDMInfo.Nfft)-1);
            knownFrequencySymbol = complex(zeros(slotOFDMInfo.Nfft,1));
            knownFrequencySymbol(slotOFDMInfo.ActiveSubcarrierIndices(dmrsSubcarriers)) = dmrsSymbolsSlots{slotIdx};
            dmrsFrequencySymbols(:,slotIdx) = ifftshift(knownFrequencySymbol);
            dmrsCyclicPrefixLengths(slotIdx) = double(slotOFDMInfo.CyclicPrefixLengths(dmrsSymbolNumber));
            if ~isempty(slotOFDMInfo.Windowing)
                dmrsWindowingSamples(slotIdx) = double(slotOFDMInfo.Windowing);
            end
        end

        if ApplyTimeFrequencyOffsets
            appliedTimingOffset = ActualTimingOffsetSamples;
            appliedFrequencyOffsetHz = ActualFrequencyOffsetHz;
        else
            appliedTimingOffset = 0;
            appliedFrequencyOffsetHz = 0;
        end

        % Guarding is expanded from the configured offset, so even offsets
        % outside the LS estimator's valid/unambiguous range do not cause
        % indexing errors. Such cases are reflected by worse BLER.
        timingCaptureGuard = max(round(waveformInfo.Nfft/2), abs(round(appliedTimingOffset)));
        [rxWaveformTB,nominalStartSample,appliedTimingOffset] = mynrApplyTimeFrequencyOffset(txWaveformTB,appliedTimingOffset,appliedFrequencyOffsetHz,waveformInfo.SampleRate,timingCaptureGuard);
        channelSampleIndex = (0:size(rxWaveformTB,1)-1).';
        noise = N0*randn(size(rxWaveformTB),"like",rxWaveformTB);
        rxWaveformTB = rxWaveformTB+noise;

        if EnableTimeFrequencyCompensation
            switch JointEstimatorDomain
                % Joint 2-D LS search over integer timing D and continuous CFO f.
                % Both domains use only known DMRS OFDM-symbol samples.
                % 二维联合LS搜索，时频偏估计
                case 'time'
                    [estimatedTimingOffset,estimatedFrequencyOffsetHz] = ...
                        mynrJointTimeDomainTimeFrequencyEstimate(rxWaveformTB,txWaveformTB,dmrsTimeSampleMask,nominalStartSample,waveformInfo.SampleRate, ...
                        timingOffsetSearchSamples,frequencyOffsetSearchRangeHz);
                % DMRS相位差分
                case 'PhaseDiff'
                    % CTO Estimation 时偏估计
                    [estimatedTimingOffset,~,timingPhaseInfo] = mynrEstimateTimingFromLSPhase(preliminaryLS,waveformInfo.Nfft,maximumPreliminaryTimingOffset);
                    % CFO Estimation 频偏估计
                    [estimatedFrequencyOffsetHz,coarseFrequencyOffsetHz,residualFrequencyOffsetHz,frequencyPhaseInfo] = ...
                        mynrEstimateFrequencyFromLSPhase(preliminaryLS,dmrsTimeSamples,waveformInfo.SampleRate,FrequencyRefinementSlotLag,frequencyOffsetSearchRangeHz);
                otherwise
                    error('JointEstimatorDomain must be timeLS or PhaseDiff.');
            end
            % Apply inverse CFO rotation to time-domain samples, then remove
            % the estimated delay. Final OFDM demodulation occurs afterward.
            frequencyCompensatedTB = rxWaveformTB .* exp(-1i*2*pi*estimatedFrequencyOffsetHz*channelSampleIndex/waveformInfo.SampleRate);
            rxWaveformCompensated = frequencyCompensatedTB(nominalStartSample+estimatedTimingOffset+(0:tbWaveformLength-1),:);
        else
            estimatedTimingOffset = 0;
            estimatedFrequencyOffsetHz = 0;
            rxWaveformCompensated = rxWaveformTB(nominalStartSample+(0:tbWaveformLength-1),:);
        end

        estimatedTimingOffsets(mc) = estimatedTimingOffset;
        estimatedFrequencyOffsetsHz(mc) = estimatedFrequencyOffsetHz;

        %% Final OFDM demodulation after timing/CFO compensation
        rxGridSlots = cell(1,NSlot);
        for nslot = 0:NSlot-1
            carrier = Carrier;
            carrier.NSlot = nslot;
            slotIdx = nslot+1;
            waveformRange = slotStartSamples(slotIdx)+ (1:slotLengths(slotIdx));
            rxGridSlots{slotIdx} = mynrOFDMDemodulate(carrier,rxWaveformCompensated(waveformRange,:));
        end

        %% Original channel estimation and decoding
        estChannelGridSlots = cell(1,NSlot);
        jointNoiseEstSlots = nVar*ones(1,NSlot);
        channelEstimateOptions = {'CDMLengths',PUSCH.DMRS.CDMLengths, 'AveragingWindow',ChannelAveragingWindow};

        if strcmpi(ChannelEstimationMode,'practical')
            % Known zero-transmit REs in DMRS symbols, determined from the
            % allocations. Non-pilot data REs must never be used as noise.
            nullIndSlots = cell(1,NSlot);
            for nullSlotIdx = 1:NSlot
                dmrsColumns = unique(floor((double(dmrsIndicesSlots{nullSlotIdx})-1)/double(Carrier.NCustomSubcarriers))+1);
                nullMask = false(Carrier.NCustomSubcarriers,Carrier.SymbolsPerSlot);
                nullMask(:,dmrsColumns) = true;
                nullMask(dmrsIndicesSlots{nullSlotIdx}) = false;
                nullMask(puschIndicesSlots{nullSlotIdx}) = false;
                nullIndSlots{nullSlotIdx} = find(nullMask);
            end
            slotGroup = 1:NSlot;
            % Full frequency/time averaging across the selected slots.
            % LS全平均信道估计
            % [averagedChannelGridGroup,~,averageChannelInfo] = ...
            %     mynrMultiSlotChannelEstimate(Carrier,rxGridSlots, dmrsIndicesSlots,dmrsSymbolsSlots,channelEstimateOptions{:});
            % 
            % % 空值噪声估计
            % % Estimate FFT-domain noise power only on known null REs.
            % [jointNoiseEst,nullNoiseInfo] = mynrNullRENoiseEstimate(rxGridSlots,nullIndSlots,dmrsIndicesSlots);
            % 
            % % 维纳插值
            % % Apply the extracted frequency/time Wiener interpolation to
            % % the fully averaged channel observations.
            % [estChannelGridGroup,wienerChannelInfo] = ...
            %     mynrMultiSlotWienerInterpolation(Carrier, averagedChannelGridGroup,dmrsIndicesSlots, ...
            %     dmrsSymbolsSlots,jointNoiseEst, 'SlotNumbers',Carrier.NSlot+slotGroup-1, ...
            %     'SlotStartSamples',nominalStartSample-1+estimatedTimingOffset+slotStartSamples);

            % Previous known-noise MMSE estimator is retained for comparison.
            % LMMSE信道估计+维纳滤波+已知噪声
            [estChannelGridGroup,jointNoiseEst,mmseChannelInfo] = ...
                mynrMultiSlotMMSEChannelEstimate(Carrier,rxGridSlots, ...
                dmrsIndicesSlots,dmrsSymbolsSlots, ...
                'NoiseVariance',nVar,'SlotNumbers',Carrier.NSlot+slotGroup-1, ...
                'SlotStartSamples',nominalStartSample-1+estimatedTimingOffset+slotStartSamples);

            % Previous raw-LS front end is retained for comparison.
            % LS信道估计+空值噪声估计+维纳滤波
            % [estChannelGridGroup,jointNoiseEst,lsChannelInfo] = ...
            %     mynrMultiSlotLSChannelEstimate(Carrier,rxGridSlots, ...
            %     dmrsIndicesSlots,dmrsSymbolsSlots, ...
            %     'NullIndSlots',nullIndSlots, ...
            %     'SlotNumbers',Carrier.NSlot+slotGroup-1, ...
            %     'SlotStartSamples',nominalStartSample-1+estimatedTimingOffset+slotStartSamples);

            estChannelGridSlots = estChannelGridGroup;
            jointNoiseEstSlots = jointNoiseEst*ones(1,NSlot);
        end

        if strcmpi(ChannelEstimationMode,'perfect')
            estChannelGridSlots = cell(1,NSlot);
            [estChannelGridSlots{:}] = deal(ones(Carrier.SubcarriersPerPRB,Carrier.SymbolsPerSlot));
            jointNoiseEstSlots = nVar*ones(1,NSlot);
        end

        ulschLLRs = zeros(GTotal,1);
        for nslot = 0:NSlot-1
            carrier = Carrier;
            carrier.NSlot = nslot;
            slotIdx = nslot+1;
            slotRange = nslot*GPerSlot+(1:GPerSlot);
            rxGrid = rxGridSlots{slotIdx};
            puschIndicesSlot = puschIndicesSlots{slotIdx};
            estChannelGrid = estChannelGridSlots{slotIdx};
            noiseEst = jointNoiseEstSlots(slotIdx);

            [puschRx,puschHest] = nrExtractResources(puschIndicesSlot,rxGrid,estChannelGrid);
            [puschEq,csi] = nrEqualizeMMSE(puschRx,puschHest,noiseEst);
            [slotLLRs,rxSymbols] = mynrPUSCHDecode(carrier,PUSCH,puschEq,noiseEst);
            csi = mynrPUSCHSymbolCSI(carrier,PUSCH,csi);
            slotLLRs = slotLLRs .* repelem(csi,nr5g.internal.getQm(PUSCH.Modulation));
            ulschLLRs(slotRange) = slotLLRs;
        end

        decodeULSCH.TransportBlockLength = trBlkSize;
        [a,blkerr(mc)] = decodeULSCH(ulschLLRs,PUSCH.Modulation,PUSCH.NumLayers,rvSeq);
    end

    BLER(2,snrIdx) = sum(blkerr)/MC;
end
toc;
time = toc;
figure;
semilogy(SNRIn,BLER(2,:));