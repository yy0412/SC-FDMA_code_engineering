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
ApplyTimeFrequencyOffsets = true;
EnableTimeFrequencyCompensation = true;
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

        %% ========== 流式信道估计/噪声估计/均衡（含 delayedInterp + 时域数据辅助） ==========
        Qm = double(nr5g.internal.getQm(PUSCH.Modulation));

        % ---- 基本量 ----
        cfg = struct();
        cfg.Nsc       = Carrier.SubcarriersPerPRB;              % 4，子载波数（本 subPRB）
        cfg.Nsym      = Carrier.SymbolsPerSlot;                 % 14，每 slot 符号数
        cfg.T_sym     = 1/(Carrier.SubcarrierSpacing*1e3);      % 符号周期（秒），120 kHz → 8.33 us
        cfg.T_slot    = cfg.Nsym * cfg.T_sym;                   % slot 周期（秒），14 符号
        cfg.deltaF    = Carrier.SubcarrierSpacing*1e3;          % 子载波间隔（Hz），120 kHz
        cfg.snrLinear = SNR;                                    % 线性 SNR，策略判据用
        cfg.Qm        = double(nr5g.internal.getQm(PUSCH.Modulation));  % 每符号比特数：1/2/4/6
        cfg.P_x_dmrs  = mean(abs(dmrsSymbolsSlots{1}).^2);      % DMRS 每 RE 功率（comb-2 增强后 = 2）
        cfg.modulation = PUSCH.Modulation;                      % 调制方式，软判决/硬判决分发用
        cfg.ctx.carrier = Carrier;                              % 供 process_one_slot 内调用 mynr* 函数
        cfg.ctx.PUSCH   = PUSCH;

        % ---- 信道估计参数 ----
        cfg.channel.model       = 'flat';       % 信道模型：'flat'（默认）| 'multipath'（需 PDP 先验）
        cfg.ch.maxWin           = 100;          % RLS 窗口物理上限（slot 数）。
        % 实际受 N_max=0.423/(|deltaFErr|·T_slot) 约束（每 slot 重算）。
        % 当前实现：N_max 截断 H_hist 长度；forgetFactor=1.0 时
        % RLS 用全部历史，未实现"硬约束"语义（暂接受）。
        cfg.ch.forgetFactor     = 1.0;          % RLS 遗忘因子。1.0 = 无遗忘（无限平均）；
        % <1 时等效窗口 N_eff=1/(1-forgetFactor)。
        % 若要实现"N_max 硬约束"，可设 forgetFactor=1-1/N_max。
        cfg.ch.warmupSlots      = 5;            % 信道 RLS 收敛判断的 slot 数，前 5 slot 标 transient
        cfg.ch.sigmaInit        = 1e4;          % RLS 初始化协方差 P=σ_init²，取大让 RLS 快速收敛
        cfg.ch.dataAided        = true;         % 是否启用数据辅助信道估计（时域增益修正）
        cfg.ch.llrThresh        = 4;            % 数据辅助信道估计的 LLR 可靠性门限，|LLR|>4 才算可靠
        cfg.ch.minAuxRE         = 8;            % 数据辅助信道估计的最小可靠 RE 数，低于则本 slot 跳过
        cfg.ch.gamma            = 0.15;         % 数据辅助信道估计的增益修正步长 [0,1]，防错误传播
        cfg.ch.holdThreshold_dB = 0.5;          % nonDMRS 策略判据：残余相位 SINR 损失 < 0.5 dB → hold
        cfg.ch.allowDelay       = true;         % 是否允许 delayedInterp（1 slot 延迟），false 时强制 hold

        % ---- 噪声估计参数 ----
        cfg.noise.gammaDMRS            = 0.3;        % DMRS 快路径 IIR 步长（样本少，取大）
        cfg.noise.gammaData            = 0.15;       % 数据慢路径 IIR 步长（样本多，取小防偏差传播）
        cfg.noise.warmupSlots          = 5;          % 噪声数据路径 converged 判断的 slot 数
        cfg.noise.llrThresh            = 4;          % 数据慢路径的 LLR 可靠性门限
        cfg.noise.minAuxRE             = 8;          % 数据慢路径的最小可靠 RE 数
        cfg.noise.inVarFloor           = 1e-6;       % inVar 下限，防负值和零
        cfg.noise.mergeMode            = 'dmrsPrimary';  % DMRS/数据合并：'dmrsOnly'|'dataOnly'|'weighted'|'dmrsPrimary'
        cfg.noise.dataCorrectionGain   = 0.3;        % dmrsPrimary 模式下数据对 DMRS 的修正增益 [0,1]
        cfg.noise.maxBiasRatio         = 0.5;        % 安全阀：偏差项/rawVar > 0.5 时回退先验
        cfg.noise.inVarPrior           = nVar;       % 首 slot 的 inVar 先验（仿真用真值；实际用长期先验）
        cfg.noise.nVarPrior            = nVar;       % 首 slot 的 nVar 先验
        cfg.noise.separationMode       = 'priorNoise';  % nVar/iVar 分离：'none'|'priorNoise'|'longTermMin'
        cfg.noise.enableDataPath       = true;       % 是否启用数据慢路径（时域残差法）
        cfg.noise.nSampleDMRS          = 2;          % DMRS 每 slot 的样本数（comb-2 共 2 点）
        cfg.noise.nSampleData          = GPerSlot/Qm;  % 数据每 slot 的样本数（52，13 数据符号×4 子载波）
        cfg.noise.betaUniformThreshold = 0.3;        % β 均匀性阈值：std(β)/mean(β)>0.3 时降权数据路径

        % ---- 时域数据辅助通用参数 ----
        cfg.dataAidedUndoCSI = true;    % LLR 已乘 csi 时，软判决前是否除回。
        % true：先除 csi 再软判决，X_soft 幅度正确；
        % false：直接软判决，X_soft 幅度被 tanh(βL) 压低。
        % 推荐 true。

        % ---- pi/2-BPSK 索引参数 ----
        cfg.pi2bpskPerSlotReset = true;   % pi/2-BPSK 的相位索引 k 是否每 slot 从 0 重置
        % （与 mynrPUSCH 每 slot 独立调制一致；当前只作语义标记）
        cfg.pi2bpskGlobalKStart = 0;      % pi/2-BPSK 全局索引起点（0-based），
        % 若发送端跨 slot 连续计数则改为 mod(nslot*NRE,2)

        % ---- 均衡参数 ----
        cfg.eq.denFloor = 1e-8;    % 均衡分母下限，防除零

        % 残余频偏（仿真注入，判据用；TB 内恒定）
        if EnableTimeFrequencyCompensation
            deltaFErr_sim = ActualFrequencyOffsetHz - estimatedFrequencyOffsetHz;
        else
            deltaFErr_sim = ActualFrequencyOffsetHz;
        end

        % ---- 状态初始化 ----
        chState = init_channel_state(cfg);
        nsState = init_noise_state(cfg);
        chState.deltaFErr     = deltaFErr_sim;
        chState.deltaFErrConf = 1.0;
        chState.fd            = 0;          % 本轮不开 Jakes

        buf = struct('rxSlotPrev',[],'dmrsIdxPrev',[],'dmrsSymPrev',[],'puschIdxPrev',[], ...
            'rxSlotCurr',[],'dmrsIdxCurr',[],'dmrsSymCurr',[],'puschIdxCurr',[]);

        ulschLLRs = zeros(GTotal,1);
        useDelayedInterpTB = false;

        % ========== 流式主循环 ==========
        for nslot = 0:NSlot-1
            slotIdx = nslot + 1;
            rxSlot  = rxGridSlots{slotIdx};
            puschIdx_slot = puschIndicesSlots{slotIdx};

            % ---- 解析 DMRS 索引与参考符号 ----
            [dmrsSub, dmrsSymNum] = ind2sub([cfg.Nsc cfg.Nsym], dmrsIndicesSlots{slotIdx});
            isDMRS3 = (dmrsSymNum == 3);
            dmrsIdx = sub2ind([cfg.Nsc cfg.Nsym], dmrsSub(isDMRS3), dmrsSymNum(isDMRS3));
            dmrsSym = dmrsSymbolsSlots{slotIdx}(isDMRS3);

            % ---- 压入 dmrsBuf（三观测环形） ----
            H_ls_n = rxSlot(dmrsIdx) ./ dmrsSym;
            chState.dmrsBuf{end+1} = H_ls_n;
            if numel(chState.dmrsBuf) > 3, chState.dmrsBuf(1) = []; end
            chState.dmrsBufCount = min(chState.dmrsBufCount + 1, 3);

            % ---- 策略判定（首 slot 时用于决定后续路径） ----
            [chState.strategy, stratInfo] = select_non_dmrs_strategy(chState, cfg);

            if nslot == 0
                % ========== 首 slot：强制 hold 输出 ==========
                [H_grid, chState, chInfo] = channel_estimate_stream_hold( ...
                    rxSlot, dmrsIdx, dmrsSym, nsState.inVar, chState, cfg);
                [slotLLRs, chState, nsState] = process_one_slot( ...
                    rxSlot, dmrsIdx, dmrsSym, puschIdx_slot, H_grid, ...
                    chState, nsState, chInfo, cfg);
                ulschLLRs(0*GPerSlot + (1:GPerSlot)) = slotLLRs;

                % 决定后续路径
                if stratInfo.sinrLossDB >= cfg.ch.holdThreshold_dB && cfg.ch.allowDelay
                    useDelayedInterpTB = true;
                end

            elseif useDelayedInterpTB
                % ========== delayedInterp 模式 ==========
                if nslot == 1
                    % dmrsBuf 不足 3，只缓存不输出
                elseif chState.dmrsBufCount >= 3
                    % 输出 slot nslot-1
                    [H_grid_nm1, chState, chInfo] = channel_interp_delayed( ...
                        chState, nsState.inVar, cfg);
                    [slotLLRs_nm1, chState, nsState] = process_one_slot( ...
                        buf.rxSlotPrev, buf.dmrsIdxPrev, buf.dmrsSymPrev, buf.puschIdxPrev, ...
                        H_grid_nm1, chState, nsState, chInfo, cfg);
                    ulschLLRs((nslot-1)*GPerSlot + (1:GPerSlot)) = slotLLRs_nm1;
                end

            else
                % ========== hold 模式 ==========
                [H_grid, chState, chInfo] = channel_estimate_stream_hold( ...
                    rxSlot, dmrsIdx, dmrsSym, nsState.inVar, chState, cfg);
                [slotLLRs, chState, nsState] = process_one_slot( ...
                    rxSlot, dmrsIdx, dmrsSym, puschIdx_slot, H_grid, ...
                    chState, nsState, chInfo, cfg);
                ulschLLRs(nslot*GPerSlot + (1:GPerSlot)) = slotLLRs;
            end

            % ---- 保存 buffer（供下一 slot delayedInterp / TB 末尾兜底用） ----
            buf.rxSlotPrev    = rxSlot;      buf.dmrsIdxPrev  = dmrsIdx;
            buf.dmrsSymPrev   = dmrsSym;     buf.puschIdxPrev = puschIdx_slot;
            buf.rxSlotCurr    = rxSlot;      buf.dmrsIdxCurr  = dmrsIdx;
            buf.dmrsSymCurr   = dmrsSym;     buf.puschIdxCurr = puschIdx_slot;
        end

        % ========== TB 末尾兜底：delayedInterp 模式下 last slot 用单 slot LS 输出 ==========
        if useDelayedInterpTB && NSlot >= 2
            H_ls_last = buf.rxSlotCurr(buf.dmrsIdxCurr) ./ buf.dmrsSymCurr;
            H_scalar  = mean(H_ls_last);                     % 标量
            H_grid_last = complex(zeros(cfg.Nsc, cfg.Nsym));
            H_grid_last(:,3)          = H_scalar;            % 4×1 ← 标量，自动扩展
            H_grid_last(:,1:2)        = H_scalar;            % 4×2 ← 标量，自动扩展  ★ 修复
            H_grid_last(:,4:cfg.Nsym) = H_scalar;            % 4×11 ← 标量，自动扩展 ★ 修复
            chInfo_last = struct('NMSE_theory', nsState.inVar / (numel(H_ls_last)*cfg.P_x_dmrs + eps), ...
                'coldStart', false, 'converged', chState.converged, ...
                'strategyUsed','hold-bootstrap', 'delaySlots', 0, ...
                'N_max', cfg.ch.maxWin, 'modelUsed','flat');
            [slotLLRs_last, chState, nsState] = process_one_slot( ...
                buf.rxSlotCurr, buf.dmrsIdxCurr, buf.dmrsSymCurr, buf.puschIdxCurr, ...
                H_grid_last, chState, nsState, chInfo_last, cfg);
            ulschLLRs((NSlot-1)*GPerSlot + (1:GPerSlot)) = slotLLRs_last;
        end

        % ========== 解码 ==========
        decodeULSCH.TransportBlockLength = trBlkSize;
        [a,blkerr(mc)] = decodeULSCH(ulschLLRs,PUSCH.Modulation,PUSCH.NumLayers,rvSeq);
    end

    BLER(2,snrIdx) = sum(blkerr)/MC;
end
toc;
time = toc;
figure;
semilogy(SNRIn,BLER(2,:));

%% ==================== 状态初始化 ====================
function s = init_channel_state(cfg)
s.H_avg      = complex(zeros(cfg.Nsc,1));
s.P          = cfg.ch.sigmaInit^2;        % flat 用标量
s.slotCount  = 0;
s.converged  = false;
s.buf_1_2    = complex(zeros(cfg.Nsc,2));
s.dmrsBuf    = {};
s.dmrsBufCount = 0;
s.H_hist     = complex(zeros(cfg.Nsc,0));
s.deltaFErr  = 0;
s.deltaFErrConf = 1.0;
s.strategy   = 'hold';
s.nmseTheory = 1.0;
s.fd         = 0;
s.pdp        = 1;
s.tau        = 0;
s.dataAidedUsed = false;
s.lastGainCorr  = 1;
end

function s = init_noise_state(cfg)
    % 噪声估计模块状态初始化
    s.inVar        = cfg.noise.inVarPrior;  % 当前总 I+N（用于下一 slot 均衡）
    s.nVar         = cfg.noise.nVarPrior;   % 热噪声分量
    s.iVar         = 0;                     % 干扰分量
    s.count        = 0;                     % data 路径的 IIR 更新计数（决定 converged）
    s.converged    = false;                 % data 路径是否 converged
    s.inVarDMRS    = cfg.noise.inVarPrior;  % DMRS 路径 IIR 状态
    s.inVarData    = cfg.noise.inVarPrior;  % 数据路径 IIR 状态
    s.inVarDMRSInit = false;                % DMRS 路径是否已初始化（首 slot 直接赋值后置 true）
    s.inVarHist    = [];                    % 最近 N 个 slot 的 inVar 历史（监控）
    s.mseXHist     = [];                    % 最近 N 个 slot 的平均判决误差（监控）
    s.biasHist     = [];                    % 最近 N 个 slot 的平均偏差（监控）
    s.betaMeanHist = [];                    % 最近 N 个 slot 的 β 均值（时域反推精度监控）
    s.dataReconUsed = false;                % 上一 slot 数据慢路径是否生效（info 输出）
end

%% ==================== 策略选择 ====================
function [strategy, info] = select_non_dmrs_strategy(state, cfg)
info = struct('degraded',false,'degradeReason','', ...
    'phiMaxDeg',0,'sinrLossDB',0);
% 残余相位影响（不依赖 buffer 状态，n=0 也可算）
T_span = 11*cfg.T_sym;
phi_max = 2*pi*state.deltaFErr*T_span;
nmse_phase = 2*(1 - sinc(phi_max/2));
sinr_loss_dB = 10*log10(1 + nmse_phase*cfg.snrLinear);
info.phiMaxDeg  = rad2deg(phi_max);
info.sinrLossDB = sinr_loss_dB;

% 首 slot 或 buffer 不足，只能 hold
if state.slotCount == 0 || state.dmrsBufCount < 3
    strategy = 'hold';
    return
end

if sinr_loss_dB < cfg.ch.holdThreshold_dB
    strategy = 'hold';
    return
end

if cfg.ch.allowDelay
    strategy = 'delayedInterp';
else
    strategy = 'hold';
    info.degraded = true;
    info.degradeReason = 'phaseErrHighNoDelay';
end
end

%% ==================== hold 路径信道估计 ====================
function [H_grid, state, info] = channel_estimate_stream_hold( ...
    rxSlot, dmrsIdx, dmrsSym, inVarPrev, state, cfg)
if state.slotCount == 0
    state.nmseTheory = 1.0;
    state.strategy   = 'hold';
    state.P          = cfg.ch.sigmaInit^2;
    state.H_avg      = complex(zeros(cfg.Nsc,1));
end
H_ls = rxSlot(dmrsIdx) ./ dmrsSym;
[H_slot_scalar, freqInfo] = channel_freq_process_flat( ...
    H_ls, inVarPrev, cfg.P_x_dmrs, cfg);
[state, H_smooth] = rls_update(H_slot_scalar, state, cfg);
state.nmseTheory = freqInfo.NMSE;

H_grid = complex(zeros(cfg.Nsc, cfg.Nsym));
H_grid(:,3)          = H_smooth;
H_grid(:,1:2)        = repmat(H_smooth, 1, 2);
H_grid(:,4:cfg.Nsym) = repmat(H_smooth, 1, cfg.Nsym - 3);

info = struct( ...
    'NMSE_theory',  freqInfo.NMSE, ...
    'coldStart',    state.slotCount == 1, ...
    'converged',    state.converged, ...
    'strategyUsed', state.strategy, ...
    'delaySlots',   0, ...
    'N_max',        min(cfg.ch.maxWin, floor(0.423/(abs(state.deltaFErr)*cfg.T_slot+eps))), ...
    'modelUsed',    'flat');
end

function [H_scalar, info] = channel_freq_process_flat(H_ls, inVar, P_x, cfg)
H_scalar = mean(H_ls);
Np = numel(H_ls);
info.NMSE = inVar / (Np * P_x + eps);
end

function [state, H_smooth] = rls_update(H_slot_scalar, state, cfg)
N_max = min(cfg.ch.maxWin, floor(0.423/(abs(state.deltaFErr)*cfg.T_slot + eps)));
if state.slotCount >= N_max
    state.P = state.P / cfg.ch.forgetFactor;
end
% flat：标量 RLS
x = 1;
y = H_slot_scalar;
K = state.P*x / (cfg.ch.forgetFactor + x'*state.P*x);
state.H_avg = state.H_avg + K*(y - x'*mean(state.H_avg));
state.P     = (state.P - K*x'*state.P) / cfg.ch.forgetFactor;
H_smooth    = state.H_avg;

state.slotCount = state.slotCount + 1;
state.converged = state.slotCount >= cfg.ch.warmupSlots;
state.H_hist    = [state.H_hist, state.H_avg];
if size(state.H_hist,2) > N_max
    state.H_hist = state.H_hist(:, end-N_max+1:end);
end
end

%% ==================== delayedInterp ====================
function [H_grid_delayed, state, info] = channel_interp_delayed(state, inVarPrev, cfg)
% state.dmrsBuf = [n-2, n-1, n]，输出 slot n-1 的 H_grid
H_nm2 = state.dmrsBuf{1};
H_nm1 = state.dmrsBuf{2};
H_n   = state.dmrsBuf{3};

[Hf_nm2, ~]  = channel_freq_process_flat(H_nm2, inVarPrev, cfg.P_x_dmrs, cfg);
[Hf_nm1, f1] = channel_freq_process_flat(H_nm1, inVarPrev, cfg.P_x_dmrs, cfg);
[Hf_n,   ~]  = channel_freq_process_flat(H_n,   inVarPrev, cfg.P_x_dmrs, cfg);

H_obs = [Hf_nm2, Hf_nm1, Hf_n];     % 1×3
t_obs = [-cfg.T_slot, 0, cfg.T_slot];

H_grid_delayed = complex(zeros(cfg.Nsc, cfg.Nsym));
for s = 1:cfg.Nsym
    if s == 3
        H_scalar = Hf_nm1;
    else
        t_s = (s - 3) * cfg.T_sym;
        w   = wiener_time_weights(t_obs, t_s, state.fd, inVarPrev, cfg);
        H_scalar = w' * H_obs.';        % 1×1
    end
    H_grid_delayed(:,s) = repmat(H_scalar, cfg.Nsc, 1);
end

% 用 n-1 的观测更新 RLS
[state, ~] = rls_update(Hf_nm1, state, cfg);
state.nmseTheory = f1.NMSE;

info = struct( ...
    'NMSE_theory',  f1.NMSE, ...
    'coldStart',    false, ...
    'converged',    state.converged, ...
    'strategyUsed','delayedInterp', ...
    'delaySlots',  1, ...
    'N_max',       min(cfg.ch.maxWin, floor(0.423/(abs(state.deltaFErr)*cfg.T_slot+eps))), ...
    'modelUsed',   'flat');
end

function w = wiener_time_weights(t_obs, t_s, fd, inVarPrev, cfg)
if fd == 0
    % 拉格朗日 3 点二次插值（无多普勒）
    w = zeros(3,1);
    w(1) = (t_s - t_obs(2))*(t_s - t_obs(3)) / ...
        ((t_obs(1)-t_obs(2))*(t_obs(1)-t_obs(3)));
    w(2) = (t_s - t_obs(1))*(t_s - t_obs(3)) / ...
        ((t_obs(2)-t_obs(1))*(t_obs(2)-t_obs(3)));
    w(3) = (t_s - t_obs(1))*(t_s - t_obs(2)) / ...
        ((t_obs(3)-t_obs(1))*(t_obs(3)-t_obs(2)));
    return;
end
R_tt = besselj(0, 2*pi*fd*abs(t_obs' - t_obs));
r_ts = besselj(0, 2*pi*fd*abs(t_obs' - t_s));
R_tt = R_tt + inVarPrev*eye(3);
w = R_tt \ r_ts;
end

%% ==================== 单 slot 完整流程封装 ====================
function [slotLLRs, chState, nsState] = process_one_slot( ...
    rxSlot, dmrsIdx, dmrsSym, puschIdx, H_grid, chState, nsState, chInfo, cfg)

carrier = cfg.ctx.carrier;
PUSCH   = cfg.ctx.PUSCH;
Qm      = cfg.Qm;

% ---- DMRS 噪声快路径 ----
[inVar_dmrs_inst, ~] = dmrs_reconstruction( ...
    rxSlot(dmrsIdx), dmrsSym, H_grid(dmrsIdx), chInfo.NMSE_theory, cfg);
nsState = iir_update(nsState, 'dmrs', inVar_dmrs_inst, cfg);

% ---- 均衡（用上一 slot 的 inVar，1 slot 延迟） ----
inVar_use = nsState.inVar;
[eqSym, csi, eqInfo] = equalize_stream( ...
    rxSlot, H_grid, inVar_use, chInfo.NMSE_theory, cfg);

% ---- 解调（含 DFT 解预编码） ----
%noise_demod = inVar_use + eqInfo.lambda;
noise_demod = inVar_use;
[slotLLRs, rxSymbols] = mynrPUSCHDecode( ...
    carrier, PUSCH, eqSym(puschIdx), noise_demod);
csi_sym = mynrPUSCHSymbolCSI(carrier, PUSCH, csi(puschIdx));
%slotLLRs = slotLLRs .* repelem(csi_sym, Qm);

% ---- 噪声慢路径（时域残差，converged 且 enableDataPath） ----
if nsState.converged && chState.converged && cfg.noise.enableDataPath
    H_data_pow = abs(H_grid(puschIdx)).^2;
    [inVar_data_inst, ~] = data_reconstruction_time_domain( ...
        rxSymbols, slotLLRs, csi_sym, H_data_pow, eqInfo.lambda, cfg);
    nsState = iir_update(nsState, 'data', inVar_data_inst, cfg);
    nsState.dataReconUsed = ~isnan(inVar_data_inst);
end
[nsState.inVar, ~] = merge_dmrs_data(nsState, cfg);
nsState.inVar = max(nsState.inVar, cfg.noise.inVarFloor);
[nsState.nVar, nsState.iVar] = separate_nVar_iVar( ...
    nsState.inVar, nsState, cfg);

% ---- 数据辅助信道估计（时域增益修正，可选） ----
if chState.converged && cfg.ch.dataAided
    [chState, ~] = data_aided_update_time( ...
        rxSymbols, slotLLRs, csi_sym, chState, cfg);
end
end

%% ==================== 噪声估计 ====================
function [inVar_inst, info] = dmrs_reconstruction(rx_dmrs, dmrsSym, H_dmrs, nmseTheory, cfg)
e = rx_dmrs - H_dmrs .* dmrsSym;
rawVar = mean(abs(e).^2);
P_x = mean(abs(dmrsSym).^2);
P_H = mean(abs(H_dmrs).^2);
bias_H = nmseTheory * P_H * P_x;
inVar_inst = max(rawVar - bias_H, cfg.noise.inVarFloor);
info = struct('fallback','','rawVar',rawVar,'biasH',bias_H, ...
    'nSample',numel(rx_dmrs));
if bias_H > cfg.noise.maxBiasRatio * rawVar
    inVar_inst = cfg.noise.inVarPrior;
    info.fallback = 'biasTooLarge';
end
end

function [inVar_inst, info] = data_reconstruction_time_domain( ...
    rxSymbols, slotLLRs, csi_sym, H_data_pow, lambda, cfg)
inVar_inst = NaN;
info = struct('used',false,'fallback','','nSample',0, ...
    'rawVar',NaN,'biasSig',NaN,'biasX',NaN, ...
    'betaMean',NaN,'P_H',NaN,'noisePow',NaN);

% 0) LLR 若已乘 csi，先除回
if cfg.dataAidedUndoCSI
    LLR_raw = slotLLRs ./ repelem(max(csi_sym, 1e-6), cfg.Qm);
else
    LLR_raw = slotLLRs;
end

% 1) 软判决（时域，pi/2-BPSK 需要 symIdx0Based）
NRE = numel(rxSymbols);
if isfield(cfg,'pi2bpskGlobalKStart')
    symIdx0Based = (0:NRE-1).' + cfg.pi2bpskGlobalKStart;
else
    symIdx0Based = (0:NRE-1).';
end
[X_soft, varX] = soft_decision_from_llr_time( ...
    LLR_raw, cfg.Qm, cfg.modulation, symIdx0Based);
X_soft = X_soft(:); varX = varX(:);
rxSymbols = rxSymbols(:); csi_sym = csi_sym(:); H_data_pow = H_data_pow(:);

if numel(X_soft) ~= numel(rxSymbols)
    error('data_reconstruction_time_domain: X_soft/rxSymbols 长度不一致 (%d vs %d)', ...
        numel(X_soft), numel(rxSymbols));
end

% 2) 高可靠筛选
llrMag = min(abs(reshape(LLR_raw, [], cfg.Qm)), [], 2);
reliable = llrMag > cfg.noise.llrThresh;
if sum(reliable) < cfg.noise.minAuxRE
    info.fallback = 'tooFewRE';
    return
end

y_r  = rxSymbols(reliable);
X_r  = X_soft(reliable);
vX_r = varX(reliable);
b_r  = csi_sym(reliable);
Hp_r = H_data_pow(reliable);

% 3) 时域残差方差
e      = y_r - X_r;
rawVar = mean(abs(e).^2);

% 4) 偏差扣除（P_X=1）
biasSig  = mean((1 - b_r).^2);
biasX    = mean(vX_r);
noisePow = max(rawVar - biasSig - biasX, cfg.noise.inVarFloor);

% 5) 时域 → 频域反推
betaMean = mean(b_r);
P_H      = mean(Hp_r);
if betaMean < 0.05
    info.fallback = 'betaTooSmall';
    return
end
sigma2_plus_lambda = noisePow * P_H / (betaMean^2);
inVar_est = max(sigma2_plus_lambda - lambda, cfg.noise.inVarFloor);

% 6) 安全阀
if ~isfinite(inVar_est) || inVar_est > cfg.noise.maxBiasRatio * 1e3 * cfg.noise.inVarPrior
    inVar_inst = cfg.noise.inVarPrior;
    info.fallback = 'biasTooLarge';
else
    inVar_inst = inVar_est;
end

info.used     = true;
info.nSample  = sum(reliable);
info.rawVar   = rawVar;
info.biasSig  = biasSig;
info.biasX    = biasX;
info.betaMean = betaMean;
info.P_H      = P_H;
info.noisePow = noisePow;
end

function [X_soft, varX] = soft_decision_from_llr_time(LLR_vec, Qm, modulation, symIdx0Based)
% 时域软判决：E[s|L] 与 Var[|s|²|L]，假设各 bit LLR 独立
% LLR 约定：L = ln(P(b=0)/P(b=1))，E[1-2b|L] = tanh(L/2)
% 输入
%   LLR_vec      : NRE×Qm 或 (NRE·Qm)×1
%   Qm           : 2, 4, 6 (BPSK 为 1)
%   modulation   : 'pi/2-BPSK' | 'QPSK' | '16QAM' | '64QAM'
%   symIdx0Based : NRE×1，时域符号 0-based 索引（pi/2-BPSK 用；其他调制可忽略）
% 输出
%   X_soft : NRE×1 复，E[s|L]
%   varX   : NRE×1 实，Var[Re(s)]+Var[Im(s)] = E[|s|²] - |E[s]|²

if nargin < 4, symIdx0Based = []; end

if isvector(LLR_vec)
    L = reshape(LLR_vec, Qm, []).';
else
    L = reshape(LLR_vec, [], Qm);
end
NRE = size(L,1);
t = tanh(L/2);     % E[1-2b_i|L]

switch modulation
    case 'pi/2-BPSK'
        if isempty(symIdx0Based)
            symIdx0Based = (0:NRE-1).';
        else
            symIdx0Based = symIdx0Based(:);
            assert(numel(symIdx0Based)==NRE, 'pi/2-BPSK: symIdx0Based 长度不匹配');
        end
        rot   = exp(1i*pi/2*mod(symIdx0Based, 2));    % ∈ {1, +j}
        X_soft = rot .* (t(:,1) / sqrt(2));
        varX   = 0.5 - 0.5*t(:,1).^2;

    case 'QPSK'
        X_soft = t(:,1)/sqrt(2) + 1i*t(:,2)/sqrt(2);
        varX   = 1 - abs(X_soft).^2;

    case '16QAM'
        t0=t(:,1); t1=t(:,2); t2=t(:,3); t3=t(:,4);
        I = t0 .* (2 - t2) / sqrt(10);
        Q = t1 .* (2 - t3) / sqrt(10);
        X_soft = I + 1i*Q;
        E_I2 = (5 - 4*t2) / 10;
        E_Q2 = (5 - 4*t3) / 10;
        varX = (E_I2 + E_Q2) - abs(X_soft).^2;

    case '64QAM'
        t0=t(:,1); t1=t(:,2); t2=t(:,3);
        t3=t(:,4); t4=t(:,5); t5=t(:,6);
        I = t0 .* (4 - t2.*(2 - t4)) / sqrt(42);
        Q = t1 .* (4 - t3.*(2 - t5)) / sqrt(42);
        X_soft = I + 1i*Q;
        E_I2 = (21 - 8*t2.*(2-t4) - 4*t4) / 42;
        E_Q2 = (21 - 8*t3.*(2-t5) - 4*t5) / 42;
        varX = (E_I2 + E_Q2) - abs(X_soft).^2;

    otherwise
        error('soft_decision_from_llr_time: 未支持调制 %s', modulation);
end
X_soft = X_soft(:);
varX   = varX(:);
end

function X_hat = hard_symbol_from_llr_time(LLR_vec, Qm, modulation, symIdx0Based)
% 从 LLR 直接硬判决重建星座点
% 约定：LLR<0 → bit=1，LLR≥0 → bit=0（与 nrSymbolDemodulate 一致）
if nargin < 4, symIdx0Based = []; end

if isvector(LLR_vec)
    L = reshape(LLR_vec, Qm, []).';
else
    L = reshape(LLR_vec, [], Qm);
end
NRE = size(L,1);
bits = L < 0;

switch modulation
    case 'pi/2-BPSK'
        if isempty(symIdx0Based)
            symIdx0Based = (0:NRE-1).';
        else
            symIdx0Based = symIdx0Based(:);
        end
        rot   = exp(1i*pi/2*mod(symIdx0Based, 2));
        base  = (1 - 2*bits(:,1)) / sqrt(2);
        X_hat = rot .* base;

    case 'QPSK'
        I = (1 - 2*bits(:,1)) / sqrt(2);
        Q = (1 - 2*bits(:,2)) / sqrt(2);
        X_hat = I + 1i*Q;

    case '16QAM'
        b0=bits(:,1); b1=bits(:,2); b2=bits(:,3); b3=bits(:,4);
        I = (1-2*b0).*(2 - (1-2*b2)) / sqrt(10);
        Q = (1-2*b1).*(2 - (1-2*b3)) / sqrt(10);
        X_hat = I + 1i*Q;

    case '64QAM'
        b0=bits(:,1); b1=bits(:,2); b2=bits(:,3);
        b3=bits(:,4); b4=bits(:,5); b5=bits(:,6);
        I = (1-2*b0).*(4 - (1-2*b2).*(2 - (1-2*b4))) / sqrt(42);
        Q = (1-2*b1).*(4 - (1-2*b3).*(2 - (1-2*b5))) / sqrt(42);
        X_hat = I + 1i*Q;

    otherwise
        error('hard_symbol_from_llr_time: 未支持调制 %s', modulation);
end
X_hat = X_hat(:);
end

function state = iir_update(state, path, inVar_inst, cfg)
% IIR 更新噪声状态
% path='dmrs'：DMRS 快路径，首 slot 直接赋值（用先验），之后 IIR 平滑
% path='data'：数据慢路径，首 slot 直接赋值，之后 IIR 平滑，并递增 count 判断 converged

    if isnan(inVar_inst), return; end

    switch path
        case 'dmrs'
            g = cfg.noise.gammaDMRS;
            if ~state.inVarDMRSInit
                % 首 slot：用本次观测直接覆盖先验（不 IIR，让估计快速起步）
                state.inVarDMRS     = inVar_inst;
                state.inVarDMRSInit = true;
            else
                % 之后：IIR 平滑
                state.inVarDMRS = (1-g)*state.inVarDMRS + g*inVar_inst;
            end
            % 注意：dmrs 路径不递增 state.count，避免干扰 data 路径的 converged 判断

        case 'data'
            g = cfg.noise.gammaData;
            if state.count == 0
                % data 路径首 slot：直接赋值
                state.inVarData = inVar_inst;
            else
                % 之后：IIR 平滑
                state.inVarData = (1-g)*state.inVarData + g*inVar_inst;
            end
            state.count = state.count + 1;
            state.converged = state.count >= cfg.noise.warmupSlots;
    end
end

function [inVar, info] = merge_dmrs_data(state, cfg)
info = struct('bias',0);
switch cfg.noise.mergeMode
    case 'dmrsOnly',    inVar = state.inVarDMRS;
    case 'dataOnly',    inVar = state.inVarData;
    case 'weighted'
        nD=cfg.noise.nSampleDMRS; nA=cfg.noise.nSampleData;
        wD=nD/(nD+nA); wA=nA/(nD+nA);
        inVar = wD*state.inVarDMRS + wA*state.inVarData;
        info.bias = abs(state.inVarDMRS - state.inVarData);
    case 'dmrsPrimary'
        delta = cfg.noise.dataCorrectionGain * (state.inVarData - state.inVarDMRS);
        inVar = state.inVarDMRS + delta;
        info.bias = abs(delta);
end
end

function [nVar, iVar] = separate_nVar_iVar(inVar, state, cfg)
switch cfg.noise.separationMode
    case 'none',        nVar = NaN; iVar = inVar;
    case 'priorNoise',  nVar = cfg.noise.nVarPrior; iVar = max(inVar - nVar, 0);
    case 'longTermMin', nVar = min([state.inVarHist, inVar]); iVar = max(inVar - nVar, 0);
end
end

%% ==================== 均衡 ====================
function [eqSym, csi, info] = equalize_stream(rx, H_grid, inVar, nmseTheory, cfg)
lambda = nmseTheory * mean(abs(H_grid(:)).^2);
den = abs(H_grid).^2 + inVar + lambda;
denClipped = den < cfg.eq.denFloor;
den(denClipped) = cfg.eq.denFloor;
eqSym = conj(H_grid).*rx./den;
csi   = min(max(abs(H_grid).^2./den, 0), 1);
bad = ~isfinite(eqSym) | ~isfinite(csi);
eqSym(bad) = 0; csi(bad) = 0;
SINReff = abs(H_grid).^2./(inVar + lambda);
info = struct('denClipped',sum(denClipped(:)),'nonFinite',sum(bad(:)), ...
    'lambda',lambda,'csiMean',mean(csi(:)),'csiStd',std(csi(:)), ...
    'SINReff_mean',mean(SINReff(:)),'SINReff_p5',prctile(SINReff(:),5));
end

%% ==================== 数据辅助信道估计（时域增益修正） ====================
function [state, info] = data_aided_update_time(rxSymbols, slotLLRs, csi_sym, state, cfg)
state.lastGainCorr  = 1;
state.dataAidedUsed = false;
info = struct('dataAidedUsed',false,'gainCorr',1, ...
    'betaTimeMean',NaN,'betaRefMean',NaN,'nAux',0);

if ~state.converged || ~cfg.ch.dataAided
    return
end

% 0) LLR 若已乘 csi，先除回
if cfg.dataAidedUndoCSI
    LLR_raw = slotLLRs ./ repelem(max(csi_sym, 1e-6), cfg.Qm);
else
    LLR_raw = slotLLRs;
end

% 1) 筛选
llrMag = min(abs(reshape(LLR_raw, [], cfg.Qm)), [], 2);
reliable = llrMag > cfg.ch.llrThresh;
if sum(reliable) < cfg.ch.minAuxRE
    return
end

% 2) 硬判决（支持所有调制；pi/2-BPSK 需要 symIdx0Based）
NRE = numel(rxSymbols);
if isfield(cfg,'pi2bpskGlobalKStart')
    symIdx0Based = (0:NRE-1).' + cfg.pi2bpskGlobalKStart;
else
    symIdx0Based = (0:NRE-1).';
end
X_hat_all = hard_symbol_from_llr_time(LLR_raw, cfg.Qm, cfg.modulation, symIdx0Based);
X_hat = X_hat_all(reliable);

% 3) 时域 β_time 与参考 β_ref
beta_time = mean(rxSymbols(reliable) ./ X_hat);
beta_ref  = mean(csi_sym(reliable));

if abs(beta_ref) < 1e-6 || ~isfinite(beta_time)
    return
end

% 4) 增益修正
gain_correction = beta_time / beta_ref;
if abs(log10(abs(gain_correction))) > 0.1
    return
end
state.H_avg = (1 - cfg.ch.gamma) * state.H_avg + ...
    cfg.ch.gamma * state.H_avg * gain_correction;

state.dataAidedUsed = true;
state.lastGainCorr  = gain_correction;

info.dataAidedUsed = true;
info.gainCorr      = gain_correction;
info.betaTimeMean  = beta_time;
info.betaRefMean   = beta_ref;
info.nAux          = sum(reliable);
end