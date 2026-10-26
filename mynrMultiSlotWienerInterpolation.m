function [hEstGrids,info] = mynrMultiSlotWienerInterpolation( ...
        carrier,channelObservationGrids,refIndSlots,refSymSlots,nVar,varargin)
%mynrMultiSlotWienerInterpolation Two-stage multi-slot Wiener interpolation.
%   [HESTGRIDS,INFO] = mynrMultiSlotWienerInterpolation(CARRIER,HOBS,
%   REFIND,REFSYM,NVAR) Wiener-interpolates channel observations at the
%   reference REs. HOBS is a cell array of channel-estimate grids, for
%   example the fully averaged output of mynrMultiSlotChannelEstimate.
%   NVAR is the received-grid noise power before division by the pilots.
%
%   Frequency interpolation is performed independently on every DMRS OFDM
%   symbol, followed by time interpolation on every occupied subcarrier.
%   Values at the reference REs are restored to their input observations;
%   Wiener estimates therefore supply only the non-reference REs.
%
%   The covariance model describes the residual timing/CFO response for
%   the configured SISO H=1 link and intentionally omits ICI.

parser = inputParser;
parser.FunctionName = 'mynrMultiSlotWienerInterpolation';
addParameter(parser,'ResidualMean',[0;0],@isnumeric);
addParameter(parser,'ResidualCovariance',[],@isnumeric);
addParameter(parser,'SlotNumbers',[],@isnumeric);
addParameter(parser,'SlotStartSamples',[],@isnumeric);
addParameter(parser,'CyclicPrefixFraction',0.5,@isnumeric);
addParameter(parser,'UseCache',true,@(x) islogical(x) && isscalar(x));
parse(parser,varargin{:});
opts = parser.Results;

validateattributes(nVar,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'}, ...
    'mynrMultiSlotWienerInterpolation','nVar');
validateattributes(opts.CyclicPrefixFraction,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',1}, ...
    'mynrMultiSlotWienerInterpolation','CyclicPrefixFraction');
if ~iscell(channelObservationGrids) || isempty(channelObservationGrids) || ...
        ~iscell(refIndSlots) || ~iscell(refSymSlots) || ...
        numel(refIndSlots) ~= numel(channelObservationGrids) || ...
        numel(refSymSlots) ~= numel(channelObservationGrids)
    error('mynrMultiSlotWienerInterpolation:SlotInputs', ...
        ['Channel-observation grids, reference indices, and symbols must ' ...
        'be matching nonempty cell arrays.']);
end

S = numel(channelObservationGrids);
K = double(carrier.NCustomSubcarriers);
L = double(carrier.SymbolsPerSlot);
if isempty(opts.SlotNumbers)
    slotNumbers = double(carrier.NSlot)+(0:S-1);
else
    slotNumbers = double(opts.SlotNumbers(:).');
end
validateattributes(slotNumbers,{'numeric'}, ...
    {'real','integer','nonnegative','finite','numel',S}, ...
    'mynrMultiSlotWienerInterpolation','SlotNumbers');
slotStarts = double(opts.SlotStartSamples(:).');
if ~isempty(slotStarts)
    validateattributes(slotStarts,{'numeric'}, ...
        {'real','integer','finite','numel',S}, ...
        'mynrMultiSlotWienerInterpolation','SlotStartSamples');
end

symbolTimes = zeros(L*S,1);
windowOffsets = zeros(L*S,1);
pilotIndex = zeros(0,1);
pilots = complex(zeros(0,1));
pilotObservations = complex(zeros(0,1));
nextSlotStart = 0;
for s = 1:S
    slotCarrier = carrier;
    slotCarrier.NSlot = slotNumbers(s);
    ofdmInfo = mynrOFDMInfo(slotCarrier);
    if s == 1
        N = double(ofdmInfo.Nfft);
        Fs = double(ofdmInfo.SampleRate);
        % Signed physical FFT bins, not resource-grid row numbers.
        bins = double(ofdmInfo.ActiveSubcarrierIndices(:))-floor(N/2)-1;
    elseif double(ofdmInfo.Nfft) ~= N || double(ofdmInfo.SampleRate) ~= Fs
        error('mynrMultiSlotWienerInterpolation:Numerology', ...
            'All input slots must use the same FFT length and sample rate.');
    end
    if isempty(opts.SlotStartSamples)
        slotStarts(s) = nextSlotStart;
    end
    nextSlotStart = slotStarts(s)+sum(double(ofdmInfo.SymbolLengths));
    cp = double(ofdmInfo.CyclicPrefixLengths(:));
    cpOffset = fix(cp*double(opts.CyclicPrefixFraction));
    startsInSlot = [0;cumsum(double(ofdmInfo.SymbolLengths(1:end-1)))];
    columns = (s-1)*L+(1:L);
    symbolTimes(columns) = slotStarts(s)+startsInSlot+cpOffset;
    windowOffsets(columns) = cpOffset-cp;

    grid = channelObservationGrids{s};
    if ~isnumeric(grid) || size(grid,1) ~= K || size(grid,2) ~= L || ...
            numel(grid) ~= K*L
        error('mynrMultiSlotWienerInterpolation:GridSizeMismatch', ...
            'Every channel-observation grid must be a K-by-L SISO matrix.');
    end
    ind = double(refIndSlots{s}(:));
    sym = double(refSymSlots{s}(:));
    validateattributes(ind,{'numeric'}, ...
        {'real','finite','integer','>=',1,'<=',K*L}, ...
        'mynrMultiSlotWienerInterpolation','refIndSlots');
    if numel(ind) ~= numel(sym) || numel(unique(ind)) ~= numel(ind)
        error('mynrMultiSlotWienerInterpolation:PilotIndices', ...
            'Each slot needs unique SISO pilot indices and one symbol per index.');
    end
    if any(~isfinite(sym)) || any(abs(sym) == 0)
        error('mynrMultiSlotWienerInterpolation:PilotSymbols', ...
            'Reference symbols must be finite and nonzero.');
    end
    observations = double(grid(ind));
    if any(~isfinite(observations))
        error('mynrMultiSlotWienerInterpolation:ChannelObservations', ...
            'Channel observations at the reference REs must be finite.');
    end
    pilotIndex = [pilotIndex;ind+(s-1)*K*L]; %#ok<AGROW>
    pilots = [pilots;sym]; %#ok<AGROW>
    pilotObservations = [pilotObservations;observations]; %#ok<AGROW>
end
if isempty(pilots)
    error('mynrMultiSlotWienerInterpolation:NoPilots', ...
        'At least one reference RE is required.');
end

allBins = repmat(bins,L*S,1);
allTimes = repelem(symbolTimes,K);
pilotColumns = floor((pilotIndex-1)/K)+1;
dmrsColumns = unique(pilotColumns,'sorted');
frequencyIndex = reshape((1:K).'+K*(dmrsColumns.'-1),[],1);
pilotPower = abs(pilots).^2;
noiseAfterDivision = double(nVar)./pilotPower;

muE = double(opts.ResidualMean(:));
validateattributes(muE,{'numeric'}, ...
    {'real','finite','numel',2}, ...
    'mynrMultiSlotWienerInterpolation','ResidualMean');
if isempty(opts.ResidualCovariance)
    [covE,fisherInfo] = localFisherPrior(N,Fs,double(nVar),pilots, ...
        allBins(pilotIndex),pilotColumns,dmrsColumns,symbolTimes,windowOffsets);
    priorSource = ['local-Fisher Gaussian approximation driven by ' ...
        'the supplied noise power'];
else
    covE = double(opts.ResidualCovariance);
    fisherInfo = struct;
    priorSource = 'externally supplied Gaussian residual moments';
end
localValidateCovariance(covE);
covE = (covE+covE.')/2;

key = struct('N',N,'Bins',allBins,'Times',allTimes, ...
    'Pilots',pilotIndex,'FrequencyIndex',frequencyIndex, ...
    'Noise',noiseAfterDivision,'Mean',muE,'Covariance',covE,'K',K);
persistent cacheKeys cacheDesigns
if isempty(cacheKeys)
    cacheKeys = {};
    cacheDesigns = {};
end
cacheHit = false;
if opts.UseCache
    for c = 1:numel(cacheKeys)
        if isequaln(key,cacheKeys{c})
            design = cacheDesigns{c};
            cacheHit = true;
            break;
        end
    end
end
if ~cacheHit
    design = localWienerDesign(key,pilotColumns,dmrsColumns);
    if opts.UseCache
        cacheKeys = [{key} cacheKeys];
        cacheDesigns = [{design} cacheDesigns];
        if numel(cacheKeys) > 4
            cacheKeys = cacheKeys(1:4);
            cacheDesigns = cacheDesigns(1:4);
        end
    end
end

innovation = pilotObservations-design.PilotMean;
estimate = design.Mean+design.Weights*innovation;
frequencyEstimate = design.Mean(frequencyIndex) ...
    +design.FrequencyWeights*innovation;

% Interpolation must not replace the supplied estimates at observed REs.
estimate(pilotIndex) = pilotObservations;
finalErrorVariance = design.ErrorVariance;
finalErrorVariance(pilotIndex) = noiseAfterDivision;
hEstGrids = cell(size(channelObservationGrids));
errorGrids = cell(size(channelObservationGrids));
for s = 1:S
    rows = (s-1)*K*L+(1:K*L);
    hEstGrids{s} = reshape(estimate(rows),K,L);
    errorGrids{s} = reshape(finalErrorVariance(rows),K,L);
end

info = struct('Method','frequency Wiener followed by time Wiener', ...
    'ObservationSource','input channel estimates at reference REs', ...
    'Target','G=A*H=A for known physical H=1', ...
    'ICIIncludedInModel',false,'NoiseVariance',double(nVar), ...
    'NoiseVarianceAfterPilotDivision',noiseAfterDivision, ...
    'PriorSource',priorSource,'ResidualMean',muE, ...
    'ResidualCovariance',covE,'FisherInfo',fisherInfo, ...
    'Nfft',N,'SampleRate',Fs,'SignedSubcarrierBins',bins, ...
    'SlotNumbers',slotNumbers,'SlotStartSamples',slotStarts, ...
    'FFTStartSamples',symbolTimes,'PilotIndices',pilotIndex, ...
    'PilotSymbols',pilots,'PilotObservations',pilotObservations, ...
    'PilotOutputUsesInputObservation',true, ...
    'PriorMean',design.Mean, ...
    'PilotPriorCovariance',design.PilotCovariance, ...
    'PilotPriorSecondMoment', ...
    design.PilotCovariance+design.PilotMean*design.PilotMean', ...
    'FinalPilotEstimates',estimate(pilotIndex), ...
    'DMRSSymbolColumns',dmrsColumns, ...
    'FrequencyStageEstimates',reshape(frequencyEstimate,K,[]), ...
    'FrequencyStageErrorVariance', ...
    reshape(design.FrequencyErrorVariance,K,[]), ...
    'FrequencyStageObservationCovariance', ...
    design.FrequencyObservationCovariance, ...
    'FrequencyWienerWeights',design.FrequencyWeights, ...
    'TimeWienerWeights',{design.TimeWeights}, ...
    'Weights',design.Weights,'ErrorVarianceGrids',{errorGrids}, ...
    'CacheHit',cacheHit);
end

function [covE,info] = localFisherPrior(N,Fs,nVar,pilots,kPilot, ...
        pilotColumns,dmrsColumns,symbolTimes,windowOffsets)
gram = zeros(2);
n = (0:N-1).';
for column = dmrsColumns.'
    selected = pilotColumns == column;
    k = kPilot(selected);
    x = pilots(selected);
    basis = exp(1i*2*pi/N*((n+windowOffsets(column))*k.'))/N;
    waveform = basis*x;
    derivativeTiming = -basis*(1i*2*pi/N*k.*x);
    derivativeHz = 1i*2*pi/Fs*(symbolTimes(column)+n).*waveform;
    J = [derivativeTiming derivativeHz];
    gram = gram+real(J'*J);
end
if nVar == 0
    covHz = zeros(2);
else
    scale = sqrt(diag(gram));
    if any(scale <= 0)
        error('mynrMultiSlotWienerInterpolation:UnidentifiablePrior', ...
            'Pilot geometry cannot determine both residual moments.');
    end
    normalized = gram./(scale*scale.');
    if min(eig(normalized)) < 1e-10
        error('mynrMultiSlotWienerInterpolation:UnidentifiablePrior', ...
            'Local Fisher matrix is singular. Supply ResidualCovariance.');
    end
    scaling = diag(1./scale);
    covHz = (nVar/N)/2*scaling*(normalized\eye(2))*scaling;
end
conversion = diag([1 2*pi/Fs]);
covE = conversion*covHz*conversion;
info = struct('ParameterUnits','[samples, Hz]', ...
    'CovarianceSamplesHz',covHz,'JacobianGram',gram, ...
    'TimeNoiseVariance',nVar/N, ...
    'Assumption','local continuous-delay, unbiased LS; known phase and H=1');
end

function localValidateCovariance(P)
validateattributes(P,{'numeric'},{'real','finite','size',[2 2]}, ...
    'mynrMultiSlotWienerInterpolation','ResidualCovariance');
if norm(P-P.','fro') > 1e-10*max(norm(P,'fro'),realmin) || any(diag(P) < 0)
    error('mynrMultiSlotWienerInterpolation:ResidualCovariance', ...
        'ResidualCovariance must be real symmetric positive semidefinite.');
end
bound = sqrt(P(1,1)*P(2,2));
if abs(P(1,2)) > bound*(1+1e-10)+1e-30
    error('mynrMultiSlotWienerInterpolation:ResidualCovariance', ...
        'Timing/frequency cross covariance exceeds the PSD bound.');
end
end

function design = localWienerDesign(key,pilotColumns,dmrsColumns)
N = key.N;
K = key.K;
k = key.Bins;
t = key.Times;
p = key.Pilots;
f = key.FrequencyIndex;
mu = localMean(N,k,t,key.Mean,key.Covariance);
muP = mu(p);
CPP = localCovariance(N,k(p),t(p),k(p),t(p), ...
    muP,muP,key.Mean,key.Covariance);
CPP = (CPP+CPP')/2;
observationCovariance = CPP+diag(key.Noise);
CFP = localCovariance(N,k(f),t(f),k(p),t(p), ...
    mu(f),muP,key.Mean,key.Covariance);
pilotCount = numel(p);
frequencyWeights = complex(zeros(numel(f),pilotCount));
for d = 1:numel(dmrsColumns)
    selected = find(pilotColumns == dmrsColumns(d));
    frequencyRows = (d-1)*K+(1:K);
    frequencyWeights(frequencyRows,selected) = localRightSolve( ...
        CFP(frequencyRows,selected), ...
        observationCovariance(selected,selected));
end
frequencyObservationCovariance = ...
    frequencyWeights*observationCovariance*frequencyWeights';
frequencyObservationCovariance = ...
    (frequencyObservationCovariance+frequencyObservationCovariance')/2;
secondMoment = localSecondMoment(N,key.Mean(2),key.Covariance(2,2));
priorVariance = max(0,secondMoment-abs(mu).^2);
frequencyErrorVariance = max(0,priorVariance(f) ...
    -real(sum(frequencyWeights.*conj(CFP),2)));

weights = complex(zeros(numel(k),pilotCount));
errorVariance = zeros(size(k));
timeWeights = cell(K,1);
for row = 1:K
    target = row:K:numel(k);
    observed = row:K:numel(f);
    rowFrequencyWeights = frequencyWeights(observed,:);
    CTP = localCovariance(N,k(target),t(target),k(p),t(p), ...
        mu(target),muP,key.Mean,key.Covariance);
    targetFrequencyCovariance = CTP*rowFrequencyWeights';
    observedCovariance = ...
        frequencyObservationCovariance(observed,observed);
    timeWeights{row} = localRightSolve( ...
        targetFrequencyCovariance,observedCovariance);
    W = timeWeights{row}*rowFrequencyWeights;
    weights(target,:) = W;
    errorVariance(target) = max(0,priorVariance(target) ...
        -2*real(sum(W.*conj(CTP),2)) ...
        +real(sum((W*observationCovariance).*conj(W),2)));
end
design = struct('Mean',mu,'PilotMean',muP, ...
    'PilotCovariance',CPP,'FrequencyWeights',frequencyWeights, ...
    'FrequencyObservationCovariance',frequencyObservationCovariance, ...
    'FrequencyErrorVariance',frequencyErrorVariance, ...
    'TimeWeights',{timeWeights},'Weights',weights, ...
    'ErrorVariance',errorVariance);
end

function mu = localMean(N,k,t,m,P)
mu = complex(zeros(size(k)));
d = -2*pi/N*k;
for n = 0:N-1
    u = t+n;
    variance = max(0,d.^2*P(1,1)+u.^2*P(2,2)+2*d.*u*P(1,2));
    mu = mu+exp(1i*(d*m(1)+u*m(2))-variance/2)/N;
end
end

function C = localCovariance(N,k1,t1,k2,t2,mu1,mu2,m,P)
if ~any(P(:))
    C = complex(zeros(numel(k1),numel(k2)));
    return;
end
deltaD = -2*pi/N*(k1-k2.');
deltaT = t1-t2.';
second = complex(zeros(size(deltaD)));
for lag = -(N-1):N-1
    u = deltaT+lag;
    variance = max(0,deltaD.^2*P(1,1)+u.^2*P(2,2) ...
        +2*deltaD.*u*P(1,2));
    second = second+(N-abs(lag))/N^2 ...
        *exp(1i*(deltaD*m(1)+u*m(2))-variance/2);
end
C = second-mu1*mu2';
end

function second = localSecondMoment(N,muOmega,varOmega)
d = (1:N-1).';
second = (N+2*sum((N-d).*exp(-d.^2*varOmega/2) ...
    .*cos(d*muOmega)))/N^2;
end

function W = localRightSolve(C,S)
% Hermitian PSD Wiener solve, including rank-deficient covariances.
S = (S+S')/2;
if ~any(S(:))
    W = complex(zeros(size(C)));
    return;
end
[V,D] = eig(S,'vector');
D = real(D);
largest = max(D);
if min(D) < -max(1,largest)*1e-9
    error('mynrMultiSlotWienerInterpolation:NonPSDCovariance', ...
        'The model observation covariance is not positive semidefinite.');
end
keep = D > max(largest*1e-11,realmin);
W = ((C*V(:,keep))./D(keep).')*V(:,keep)';
end
