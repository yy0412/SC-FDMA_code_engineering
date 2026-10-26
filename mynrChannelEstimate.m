function [hEst,nVar,info] = mynrChannelEstimate(carrier,rxGrid,refInd,refSym,varargin)
%mynrChannelEstimate Practical channel estimation on a custom-K grid.
%
%   The estimator follows the same high-level flow as nrChannelEstimate:
%   LS estimates on reference REs, CDM despreading/averaging, frequency
%   interpolation, noise estimation, and time interpolation. All operations
%   are performed directly on K = NSizeGrid*SubcarriersPerPRB subcarriers.

[K,N,R] = size(rxGrid);
if K ~= carrier.NCustomSubcarriers
    error('mynrChannelEstimate:GridSizeMismatch', ...
        'rxGrid first dimension must match carrier.NCustomSubcarriers.');
end

opts = localParseOptions(varargin{:}); % Keep the orginal toolbox optional input check
[refInd,refSym,P] = localPrepareReferenceInputs(K,N,refInd,refSym);  % Keep the orginal toolbox optional input check
ofdmInfo = mynrOFDMInfo(carrier);
cirHalfWindow = max(1,floor(min(ofdmInfo.CyclicPrefixLengths)/ofdmInfo.Nfft*K));

hRef = complex(zeros(K,N,R,P,'like',rxGrid));
hRef(:) = NaN;
noiseResiduals = complex(zeros(0,1,'like',rxGrid));

% ---------------------------------------------------------------------
% LS estimation, CDM despreading
% ---------------------------------------------------------------------
if ~isempty(refInd)
    [kAll,nAll,pAll] = ind2sub([K N P],refInd);
    
    for p = 1:P % For each port
        nSet = unique(nAll(pAll==p)).'; % The DMRS symbol locations for all slots
        for ni = 1:numel(nSet) 
            n = nSet(ni);
            this = pAll==p & nAll==n;
            k = kAll(this);
            sym = refSym(this);

            for r = 1:R
                rx = rxGrid(:,:,r);
                hLS = rx(sub2ind([K N],k,repmat(n,size(k)))) ./ sym; % LS estimation
                [hCDM,res] = localDespreadFD(hLS,opts.CDMLengths(1)); % No FD-CDM
                noiseResiduals = [noiseResiduals; res(:)]; %#ok<AGROW>
                hRef(:,n,r,p) = localAssignObserved(hRef(:,n,r,p),k,hCDM);
            end
        end
    end

    if opts.CDMLengths(2) > 1 % No TD-CDM
        [hRef,res] = localDespreadTD(hRef,opts.CDMLengths(2));
        noiseResiduals = [noiseResiduals; res(:)];
    end
end

%avgWindow = localResolveAveragingWindow(opts.AveragingWindow,K,hRef); % Provide the average window for both directions
avgWindow(1) = 2*carrier.NCustomSubcarriers;
%noiseResiduals = [noiseResiduals; localAveragingResiduals(hRef,avgWindow(1))]; % The noise elements

hEst = complex(zeros(K,N,R,P,'like',rxGrid));
for p = 1:P
    for r = 1:R
        dmrsSymbols = find(squeeze(any(isfinite(hRef(:,:,r,p)),1))).';
        if isempty(dmrsSymbols)
            hEst(:,:,r,p) = 1;
            continue;
        end

        hAtRefSymbols = complex(zeros(K,numel(dmrsSymbols),'like',rxGrid));
        for i = 1:numel(dmrsSymbols)
            n = dmrsSymbols(i);
            hAtRefSymbols(:,i) = localFrequencyLine(hRef(:,n,r,p),avgWindow(1),cirHalfWindow); % Interpolation in frequency domain
        end
        avgWindow(2) =2*numel(dmrsSymbols);
        hAtRefSymbols = localAverageTimeLines(hAtRefSymbols,avgWindow(2));% Moving average in time domain at the DMRS location
        hEst(:,:,r,p) = localInterpolateTime(hAtRefSymbols,dmrsSymbols,N);% Linear interpolation in time domain
    end
end

if strcmpi(opts.Interpolation,'off')
    hEst = localBlankNonReferenceREs(hEst,isfinite(hRef));
end
% ---------------------------------------------------------------------
% Noise estimation
% ---------------------------------------------------------------------
noiseResiduals = [noiseResiduals; localAveragingResiduals(hRef,hEst,sym,rxGrid)]; % The noise elements
nVar = localNoiseVariance(noiseResiduals);
info = struct('AveragingWindow',avgWindow, ...
    'CDMLengths',opts.CDMLengths, ...
    'CustomSubcarriers',K);
end

function opts = localParseOptions(varargin)
opts = struct('CDMLengths',[1 1], ...
    'AveragingWindow',[0 0], ...
    'Interpolation','on');

if mod(numel(varargin),2) ~= 0
    error('mynrChannelEstimate:InvalidNameValue', ...
        'Optional inputs must be supplied as name-value pairs.');
end

for i = 1:2:numel(varargin)
    name = lower(char(varargin{i}));
    value = varargin{i+1};
    switch name
        case 'cdmlengths'
            validateattributes(value,{'numeric'},{'real','integer','numel',2,'>=',1}, ...
                'mynrChannelEstimate','CDMLengths');
            opts.CDMLengths = double(reshape(value,1,2));
        case 'averagingwindow'
            validateattributes(value,{'numeric'},{'real','integer','numel',2,'>=',0}, ...
                'mynrChannelEstimate','AveragingWindow');
            opts.AveragingWindow = double(reshape(value,1,2));
        case 'interpolation'
            opts.Interpolation = validatestring(value,{'on','off'}, ...
                'mynrChannelEstimate','Interpolation');
        otherwise
            error('mynrChannelEstimate:UnknownOption', ...
                'Unknown option "%s".',varargin{i});
    end
end
end

function [refInd,refSym,P] = localPrepareReferenceInputs(K,N,refInd,refSym)
refInd = double(refInd(:));
refSym = refSym(:);

if numel(refInd) ~= numel(refSym)
    error('mynrChannelEstimate:UnequalSymIndCount', ...
        'REFIND and REFSYM must contain the same number of elements.');
end

valid = refInd > 0 & refSym ~= 0;
refInd = refInd(valid);
refSym = refSym(valid);

if isempty(refInd)
    P = 1;
else
    P = ceil(max(refInd)/(K*N));
end
end

function [hOut,res] = localDespreadFD(hIn,fdCDM)
hOut = hIn;
res = complex(zeros(0,size(hIn,2),'like',hIn));

if fdCDM <= 1 || isempty(hIn)
    return;
end

for first = 1:fdCDM:numel(hIn)
    group = first:min(first+fdCDM-1,numel(hIn));
    hAvg = mean(hIn(group,:),1);
    res = [res; hIn(group,:) - hAvg]; %#ok<AGROW>
    hOut(group,:) = repmat(hAvg,numel(group),1);
end
end

function hLine = localAssignObserved(hLine,k,h)
[ku,~,ic] = unique(k,'stable');
for i = 1:numel(ku)
    hLine(ku(i)) = mean(h(ic==i));
end
end

function [hRef,res] = localDespreadTD(hRef,tdCDM)
[K,~,R,P] = size(hRef);
res = complex(zeros(0,1,'like',hRef));

for p = 1:P
    for r = 1:R
        for k = 1:K
            n = find(isfinite(hRef(k,:,r,p)));
            if isempty(n)
                continue;
            end
            for first = 1:tdCDM:numel(n)
                group = n(first:min(first+tdCDM-1,numel(n)));
                h = hRef(k,group,r,p);
                hAvg = mean(h);
                res = [res; h(:) - hAvg]; %#ok<AGROW>
                hRef(k,group,r,p) = hAvg;
            end
        end
    end
end
end

function avgWindow = localResolveAveragingWindow(requested,K,hRef)
avgWindow = requested;
if avgWindow(1) == 0
    avgWindow(1) = min(localOddFloor(K),23);
end
if avgWindow(2) == 0
    nRefSymbols = 0;
    for n = 1:size(hRef,2)
        if any(isfinite(hRef(:,n,:,:)),'all')
            nRefSymbols = nRefSymbols + 1;
        end
    end
    if nRefSymbols <= 1
        avgWindow(2) = 1;
    else
        avgWindow(2) = min(localOddFloor(nRefSymbols),3);
    end
end
avgWindow = max(avgWindow,1);
avgWindow(1) = localForceOdd(avgWindow(1));
avgWindow(2) = localForceOdd(avgWindow(2));
end

function n = localOddFloor(n)
n = max(1,2*floor((double(n)-1)/2)+1);
end

function n = localForceOdd(n)
n = double(n);
if mod(n,2) == 0
    n = max(1,n-1);
end
end

% function res = localAveragingResiduals(hRef,freqWindow)
% [~,N,R,P] = size(hRef);
% res = complex(zeros(0,1,'like',hRef));
% for p = 1:P
%     for r = 1:R
%         for n = 1:N
%             k = find(isfinite(hRef(:,n,r,p)));
%             if numel(k) < 2 || freqWindow <= 1
%                 continue;
%             end
%             h = hRef(k,n,r,p);
%             hAvg = localMovingAverage(h,freqWindow); % Moving average on frequency direction
%             res = [res; h(:) - hAvg(:)]; %#ok<AGROW>
%         end
%     end
% end
% end

function res = localAveragingResiduals(hRef,hEst,sym,rxGrid)
[~,N,R,P] = size(hRef);
res = complex(zeros(0,1,'like',hRef));
for p = 1:P
    for r = 1:R
        for n = 1:N
            k = find(isfinite(hRef(:,n,r,p)));
            if numel(k) < 2 
                continue;
            end
            hAvg = hEst(k,n,r,p); % Extract the hEst on the DMRS location 
            res = [res; rxGrid(k,n,r,p) - hAvg.*sym]; % y = h*x + n ---> n = y - hEst*x
        end
    end
end
end

function hLine = localFrequencyLine(hObserved,freqWindow,cirHalfWindow)
K = numel(hObserved);
k = find(isfinite(hObserved));
if isempty(k)
    hLine = ones(K,1,'like',hObserved);
    return;
end

h = hObserved(k);
if freqWindow > 1
    h = localMovingAverage(h,freqWindow); % Moving average in frequency domain
end

if isscalar(k)
    hLine = repmat(h,K,1);
else
    [k,ia] = unique(k,'stable');
    h = h(ia);
    hLine = interp1(k,h,(1:K).','linear','extrap');
end
hLine = localDenoiseCIR(hLine,cirHalfWindow);
end

function hLine = localDenoiseCIR(hLine,cirHalfWindow)
K = numel(hLine);
if K <= 2
    return;
end

cirHalfWindow = min(double(cirHalfWindow),floor((K-1)/2));
mask = false(K,1);
mask(1:cirHalfWindow+1) = true;
mask(end-cirHalfWindow+1:end) = true;

cir = ifft(ifftshift(hLine),K);
hLine = fftshift(fft(cir.*mask,K));
end

function hAvg = localMovingAverage(h,window)
hAvg = h;
halfWindow = floor(double(window)/2);
for i = 1:numel(h)
    first = max(1,i-halfWindow); % The upper-half window
    last = min(numel(h),i+halfWindow); % The lower-half window
    hAvg(i) = mean(h(first:last)); % The average
end
end

function hLines = localAverageTimeLines(hLines,timeWindow)
% Moving average in time domain
if timeWindow <= 1 || size(hLines,2) <= 1
    return;
end

halfWindow = floor(double(timeWindow)/2);
hIn = hLines;
for i = 1:size(hLines,2)
    first = max(1,i-halfWindow);
    last = min(size(hLines,2),i+halfWindow);
    hLines(:,i) = mean(hIn(:,first:last),2);
end
end

function h = localInterpolateTime(hAtRefSymbols,dmrsSymbols,N)
K = size(hAtRefSymbols,1);
if isscalar(dmrsSymbols)
    h = repmat(hAtRefSymbols,1,N);
else
    h = complex(zeros(K,N,'like',hAtRefSymbols));
    for k = 1:K
        h(k,:) = interp1(dmrsSymbols,hAtRefSymbols(k,:),1:N,'linear','extrap');
    end
end
end

function hEst = localBlankNonReferenceREs(hEst,refMask)
hEst(~refMask) = 0;
end

function nVar = localNoiseVariance(res)
res = res(isfinite(real(res)) & isfinite(imag(res)));
if isempty(res)
    nVar = eps;
else 
    nVar = max(var(res),eps);
end
end
