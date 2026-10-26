function [hEstGrids,nVar,info] = mynrMultiSlotMMSEChannelEstimate( ...
        carrier,rxGridSlots,refIndSlots,refSymSlots,varargin)
%mynrMultiSlotMMSEChannelEstimate SISO H=1, residual-sync LMMSE/Wiener.
%   Estimates G(p,k)=A(p,k)*H(k)=A(p,k), omitting ICI in the observation
%   model R=X*G+W. All statistics use the complete N-sample Dirichlet factor:
%     A(p,k)=mean_n exp(j*(-2*pi*k/N*eD+(T(p)+n)*eOmega)).
%   e=[eD;eOmega] is modelled as a real jointly Gaussian residual vector.
%   This is a local approximation for an integer timing-search receiver.
%
%   Required name/value:
%     NoiseVariance       FFT-domain E|W|^2 BEFORE division by the pilots.
%   Optional name/value:
%     ResidualMean        [samples; rad/sample], default [0;0].
%     ResidualCovariance  2-by-2 covariance in the same units. If empty,
%                         use a local Fisher approximation from the known
%                         DMRS waveform and NoiseVariance/Nfft. This is not
%                         a measured covariance or a discrete-error law.
%     SlotNumbers         Actual carrier NSlot values for the input slots.
%     SlotStartSamples    Absolute zero-based CP starts of the COMPENSATED
%                         slots in the original capture, including Dhat.
%                         Default: concatenate slots beginning at sample 0.
%     CyclicPrefixFraction Default 0.5, matching mynrOFDMDemodulate.
%     UseCache            Cache up to four geometry/prior Wiener designs.
%
%   Frequency: for each DMRS OFDM symbol, jointly estimate the pilot REs
%   and Wiener-interpolate all occupied subcarriers from its raw R./X.
%   Time: for each subcarrier, Wiener-interpolate all OFDM symbols from
%   those frequency-stage estimates. Their FULL observation covariance
%   and target cross-covariance are propagated; they are not treated as
%   independent pilots with an added posterior-error noise variance.
%   Neither stage uses linear/spline interpolation or a Kalman filter.
%
%   H=1 is known; G is not set to 1 except for a deterministic zero-error
%   prior. The received data REs and the true injected offsets are unused.
%   NVAR is returned unchanged in received-RE units for the equalizer.
%   INFO.ErrorVarianceGrids instead describes model estimation uncertainty.
%   Same-data synchronization can correlate residual errors and noise;
%   the standard independent-noise Wiener model here is an approximation.

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser,'NoiseVariance',[],@isnumeric);
addParameter(parser,'ResidualMean',[0;0],@isnumeric);
addParameter(parser,'ResidualCovariance',[],@isnumeric);
addParameter(parser,'SlotNumbers',[],@isnumeric);
addParameter(parser,'SlotStartSamples',[],@isnumeric);
addParameter(parser,'CyclicPrefixFraction',0.5,@isnumeric);
addParameter(parser,'UseCache',true,@(x) islogical(x) && isscalar(x));
parse(parser,varargin{:});
opts = parser.Results;
validateattributes(opts.NoiseVariance,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,'NoiseVariance');
nVar = double(opts.NoiseVariance);
validateattributes(opts.CyclicPrefixFraction,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',1},mfilename,'CyclicPrefixFraction');
if ~iscell(rxGridSlots) || isempty(rxGridSlots) || ...
        ~iscell(refIndSlots) || ~iscell(refSymSlots) || ...
        numel(refIndSlots) ~= numel(rxGridSlots) || ...
        numel(refSymSlots) ~= numel(rxGridSlots)
    error('mynrMultiSlotMMSEChannelEstimate:SlotInputs', ...
        'Received grids, reference indices and symbols must be matching nonempty cells.');
end
S = numel(rxGridSlots);
K = double(carrier.NCustomSubcarriers);
L = double(carrier.SymbolsPerSlot);
if isempty(opts.SlotNumbers)
    slotNumbers = double(carrier.NSlot)+(0:S-1);
else
    slotNumbers = double(opts.SlotNumbers(:).');
end
validateattributes(slotNumbers,{'numeric'}, ...
    {'real','integer','nonnegative','finite','numel',S},mfilename,'SlotNumbers');
slotStarts = double(opts.SlotStartSamples(:).');
if ~isempty(slotStarts)
    validateattributes(slotStarts,{'numeric'}, ...
        {'real','integer','finite','numel',S},mfilename,'SlotStartSamples');
end

symbolTimes = zeros(L*S,1);
windowOffsets = zeros(L*S,1);
pilotIndex = zeros(0,1);
pilots = complex(zeros(0,1));
receivedPilots = complex(zeros(0,1));
nextSlotStart = 0;
for s = 1:S
    slotCarrier = carrier;
    slotCarrier.NSlot = slotNumbers(s);
    o = mynrOFDMInfo(slotCarrier);
    if s == 1
        N = double(o.Nfft);
        Fs = double(o.SampleRate);
        % Signed physical bins; do not use resource-grid rows as FFT bins.
        bins = double(o.ActiveSubcarrierIndices(:))-floor(N/2)-1;
    elseif double(o.Nfft) ~= N || double(o.SampleRate) ~= Fs
        error('mynrMultiSlotMMSEChannelEstimate:Numerology', ...
            'All input slots must use the same FFT length and sample rate.');
    end
    if isempty(opts.SlotStartSamples)
        slotStarts(s) = nextSlotStart;
    end
    nextSlotStart = slotStarts(s)+sum(double(o.SymbolLengths));
    cp = double(o.CyclicPrefixLengths(:));
    cpOffset = fix(cp*double(opts.CyclicPrefixFraction));
    startsInSlot = [0;cumsum(double(o.SymbolLengths(1:end-1)))];
    columns = (s-1)*L+(1:L);
    symbolTimes(columns) = slotStarts(s)+startsInSlot+cpOffset;
    windowOffsets(columns) = cpOffset-cp;

    grid = rxGridSlots{s};
    validateattributes(grid,{'numeric'},{'2d','size',[K L]},mfilename,'rxGridSlots');
    ind = double(refIndSlots{s}(:));
    sym = double(refSymSlots{s}(:));
    validateattributes(ind,{'numeric'}, ...
        {'real','finite','integer','>=',1,'<=',K*L},mfilename,'refIndSlots');
    if numel(ind) ~= numel(sym) || numel(unique(ind)) ~= numel(ind)
        error('mynrMultiSlotMMSEChannelEstimate:PilotIndices', ...
            'Each slot needs unique SISO pilot indices and one symbol per index.');
    end
    if any(~isfinite(sym)) || any(abs(sym) == 0)
        error('mynrMultiSlotMMSEChannelEstimate:PilotSymbols', ...
            'Reference symbols must be finite and nonzero.');
    end
    y = double(grid(ind));
    if any(~isfinite(y))
        error('mynrMultiSlotMMSEChannelEstimate:ReceivedPilots', ...
            'Received pilot observations must be finite.');
    end
    pilotIndex = [pilotIndex; ind+(s-1)*K*L]; %#ok<AGROW>
    pilots = [pilots; sym]; %#ok<AGROW>
    receivedPilots = [receivedPilots; y(:)]; %#ok<AGROW>
end
if isempty(pilots)
    error('mynrMultiSlotMMSEChannelEstimate:NoPilots','At least one DMRS RE is required.');
end
allBins = repmat(bins,L*S,1);
allTimes = repelem(symbolTimes,K);
pilotColumns = floor((pilotIndex-1)/K)+1;
dmrsColumns = unique(pilotColumns,'sorted');
frequencyIndex = reshape((1:K).'+K*(dmrsColumns.'-1),[],1);
muE = double(opts.ResidualMean(:));
validateattributes(muE,{'numeric'}, ...
    {'real','finite','numel',2},mfilename,'ResidualMean');
if isempty(opts.ResidualCovariance)
    [covE,fisherInfo] = localFisherPrior(N,Fs,nVar,pilots, ...
        allBins(pilotIndex),pilotColumns,dmrsColumns,symbolTimes,windowOffsets);
    priorSource = 'local-Fisher Gaussian approximation; not calibrated integer-error statistics';
else
    covE = double(opts.ResidualCovariance);
    fisherInfo = struct;
    priorSource = 'externally supplied Gaussian residual moments';
end
localValidateCovariance(covE);
covE = (covE+covE.')/2;

pilotPower = abs(pilots).^2;
z = receivedPilots.*conj(pilots)./pilotPower;
noiseAfterDivision = nVar./pilotPower;
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
    design = localDesign(key,pilotColumns,dmrsColumns);
    if opts.UseCache
        cacheKeys = [{key} cacheKeys];
        cacheDesigns = [{design} cacheDesigns];
        if numel(cacheKeys) > 4
            cacheKeys = cacheKeys(1:4);
            cacheDesigns = cacheDesigns(1:4);
        end
    end
end
innovation = z-design.PilotMean;
estimate = design.Mean+design.Weights*innovation;
frequencyEstimate = design.Mean(frequencyIndex)+design.FrequencyWeights*innovation;
hEstGrids = cell(size(rxGridSlots));
errorGrids = cell(size(rxGridSlots));
for s = 1:S
    rows = (s-1)*K*L+(1:K*L);
    hEstGrids{s} = reshape(estimate(rows),K,L);
    errorGrids{s} = reshape(design.ErrorVariance(rows),K,L);
end

info = struct('Method','frequency Wiener followed by time Wiener', ...
    'Target','G=A*H=A for known physical H=1', ...
    'ICIIncludedInModel',false,'NoiseVariance',nVar, ...
    'NoiseVarianceAfterPilotDivision',noiseAfterDivision, ...
    'PriorSource',priorSource,'ResidualMean',muE, ...
    'ResidualCovariance',covE,'FisherInfo',fisherInfo, ...
    'Nfft',N,'SampleRate',Fs,'SignedSubcarrierBins',bins, ...
    'SlotNumbers',slotNumbers,'SlotStartSamples',slotStarts, ...
    'FFTStartSamples',symbolTimes,'PilotIndices',pilotIndex, ...
    'PilotSymbols',pilots,'PilotLS',z, ...
    'PriorMean',design.Mean,'PilotPriorCovariance',design.PilotCovariance, ...
    'PilotPriorSecondMoment',design.PilotCovariance+design.PilotMean*design.PilotMean', ...
    'PilotLMMSE',design.PilotMean+design.PilotWeights*innovation, ...
    'FinalPilotEstimates',estimate(pilotIndex), ...
    'DMRSSymbolColumns',dmrsColumns, ...
    'FrequencyStageEstimates',reshape(frequencyEstimate,K,[]), ...
    'FrequencyStageErrorVariance',reshape(design.FrequencyErrorVariance,K,[]), ...
    'FrequencyStageObservationCovariance',design.FrequencyObservationCovariance, ...
    'FrequencyWienerWeights',design.FrequencyWeights, ...
    'TimeWienerWeights',{design.TimeWeights},'Weights',design.Weights, ...
    'ErrorVarianceGrids',{errorGrids},'CacheHit',cacheHit);
end

function [covE,info] = localFisherPrior(N,Fs,nVar,pilots,kPilot, ...
        pilotColumns,dmrsColumns,symbolTimes,windowOffsets)
% Jacobian of known H=1 pilot-only waveforms, evaluated at local alignment.
% Work in [sample, Hz] coordinates for conditioning, then convert to omega.
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
        error('mynrMultiSlotMMSEChannelEstimate:UnidentifiablePrior', ...
            'Pilot geometry cannot determine both residual moments. Supply ResidualCovariance.');
    end
    normalized = gram./(scale*scale.');
    if min(eig(normalized)) < 1e-10
        error('mynrMultiSlotMMSEChannelEstimate:UnidentifiablePrior', ...
            'Local Fisher matrix is singular. Supply calibrated ResidualCovariance.');
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
    mfilename,'ResidualCovariance');
if norm(P-P.','fro') > 1e-10*max(norm(P,'fro'),realmin) || any(diag(P) < 0)
    error('mynrMultiSlotMMSEChannelEstimate:ResidualCovariance', ...
        'ResidualCovariance must be real symmetric positive semidefinite.');
end
bound = sqrt(P(1,1)*P(2,2));
if abs(P(1,2)) > bound*(1+1e-10)+1e-30
    error('mynrMultiSlotMMSEChannelEstimate:ResidualCovariance', ...
        'Timing/frequency cross covariance exceeds the PSD bound.');
end
end

function design = localDesign(key,pilotColumns,dmrsColumns)
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
S = CPP+diag(key.Noise);
CFP = localCovariance(N,k(f),t(f),k(p),t(p), ...
    mu(f),muP,key.Mean,key.Covariance);
Pcount = numel(p);
B = complex(zeros(numel(f),Pcount));
pilotWeights = complex(zeros(Pcount));
for d = 1:numel(dmrsColumns)
    selected = find(pilotColumns == dmrsColumns(d));
    frequencyRows = (d-1)*K+(1:K);
    B(frequencyRows,selected) = localRightSolve(CFP(frequencyRows,selected),S(selected,selected));
    pilotWeights(selected,selected) = localRightSolve(CPP(selected,selected),S(selected,selected));
end
% fhat-muF=B*(z-muP). Its covariance is B*S*B', not CF+posterior error.
frequencyObservationCovariance = B*S*B';
frequencyObservationCovariance = ...
    (frequencyObservationCovariance+frequencyObservationCovariance')/2;
secondMoment = localSecondMoment(N,key.Mean(2),key.Covariance(2,2));
priorVariance = max(0,secondMoment-abs(mu).^2);
frequencyErrorVariance = max(0,priorVariance(f)-real(sum(B.*conj(CFP),2)));
weights = complex(zeros(numel(k),Pcount));
errorVariance = zeros(size(k));
timeWeights = cell(K,1);
for row = 1:K
    target = row:K:numel(k);
    observed = row:K:numel(f);
    Bk = B(observed,:);
    CTP = localCovariance(N,k(target),t(target),k(p),t(p), ...
        mu(target),muP,key.Mean,key.Covariance);
    CTargetObserved = CTP*Bk';
    observedCovariance = frequencyObservationCovariance(observed,observed);
    timeWeights{row} = localRightSolve(CTargetObserved,observedCovariance);
    W = timeWeights{row}*Bk;
    weights(target,:) = W;
    errorVariance(target) = max(0,priorVariance(target) ...
        -2*real(sum(W.*conj(CTP),2))+real(sum((W*S).*conj(W),2)));
end
design = struct('Mean',mu,'PilotMean',muP,'PilotCovariance',CPP, ...
    'PilotWeights',pilotWeights,'FrequencyWeights',B, ...
    'FrequencyObservationCovariance',frequencyObservationCovariance, ...
    'FrequencyErrorVariance',frequencyErrorVariance, ...
    'TimeWeights',{timeWeights},'Weights',weights,'ErrorVariance',errorVariance);
end

function mu = localMean(N,k,t,m,P)
mu = complex(zeros(size(k)));
d = -2*pi/N*k;
for n = 0:N-1
    u = t+n;
    v = max(0,d.^2*P(1,1)+u.^2*P(2,2)+2*d.*u*P(1,2));
    mu = mu+exp(1i*(d*m(1)+u*m(2))-v/2)/N;
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
% n-m=d occurs N-|d| times. This exactly reduces N^2 terms to 2N-1.
for lag = -(N-1):N-1
    u = deltaT+lag;
    v = max(0,deltaD.^2*P(1,1)+u.^2*P(2,2)+2*deltaD.*u*P(1,2));
    second = second+(N-abs(lag))/N^2*exp(1i*(deltaD*m(1)+u*m(2))-v/2);
end
C = second-mu1*mu2';
end

function second = localSecondMoment(N,muOmega,varOmega)
d = (1:N-1).';
second = (N+2*sum((N-d).*exp(-d.^2*varOmega/2).*cos(d*muOmega)))/N^2;
end

function W = localRightSolve(C,S)
% Hermitian PSD Wiener solve, including deterministic/rank-deficient priors.
S = (S+S')/2;
if ~any(S(:))
    W = complex(zeros(size(C)));
    return;
end
[V,D] = eig(S,'vector');
D = real(D);
largest = max(D);
if min(D) < -max(1,largest)*1e-9
    error('mynrMultiSlotMMSEChannelEstimate:NonPSDCovariance', ...
        'The model observation covariance is not positive semidefinite.');
end
keep = D > max(largest*1e-11,realmin);
W = ((C*V(:,keep))./D(keep).')*V(:,keep)';
end
