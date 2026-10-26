function [timingOffset,floatingPointDelay,info] = mynrEstimateTimingFromLSPhase( ...
        preliminaryLS,nfft,maxAbsTimingOffset)
%mynrEstimateTimingFromLSPhase Joint timing from same-slot DMRS phase pairs.
%   Every pilot pair contributes C=Z(k2,s)*conj(Z(k1,s)). Following the
%   phase-difference derivation, sum raw C over slots with the same
%   subcarrier separation, then take its angle. If there are several
%   separations, average their delay estimates with equal coefficients.
%   No channel-gain fit, search, or extra correlation-magnitude weights.

validateattributes(preliminaryLS,{'numeric'},{'2d'},mfilename,'preliminaryLS');
validateattributes(nfft,{'numeric'}, ...
    {'scalar','real','finite','positive','integer'},mfilename,'nfft');
validateattributes(maxAbsTimingOffset,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,'maxAbsTimingOffset');

[~,numberOfSlots] = size(preliminaryLS);
pairProducts = complex(zeros(0,1));
separations = zeros(0,1);
slotNumbers = zeros(0,1);
for s = 1:numberOfSlots
    rows = find(isfinite(preliminaryLS(:,s)) ...
        & abs(preliminaryLS(:,s)) > 0);
    for a = 1:numel(rows)-1
        for b = a+1:numel(rows)
            product = preliminaryLS(rows(b),s) ...
                *conj(preliminaryLS(rows(a),s));
            pairProducts(end+1,1) = product; %#ok<AGROW>
            separations(end+1,1) = rows(b)-rows(a); %#ok<AGROW>
            slotNumbers(end+1,1) = s; %#ok<AGROW>
        end
    end
end

if isempty(pairProducts)
    timingOffset = 0;
    floatingPointDelay = 0;
    info = struct('PairCount',0,'Method','unweighted phase difference');
    return;
end

uniqueSeparations = unique(separations);
groupProducts = complex(zeros(numel(uniqueSeparations),1));
groupDelays = zeros(numel(uniqueSeparations),1);
for g = 1:numel(uniqueSeparations)
    separation = uniqueSeparations(g);
    groupProducts(g) = sum(pairProducts(separations == separation));
    groupDelays(g) = -double(nfft)*angle(groupProducts(g)) ...
        /(2*pi*separation);
end
floatingPointDelay = mean(groupDelays);

timingOffset = round(floatingPointDelay);
timingOffset = max(-floor(double(maxAbsTimingOffset)), ...
    min(floor(double(maxAbsTimingOffset)),timingOffset));
info = struct('Method','unweighted phase difference', ...
    'PairCount',numel(pairProducts),'PairSlotNumbers',slotNumbers, ...
    'SubcarrierSeparations',separations, ...
    'PairProducts',pairProducts, ...
    'UniqueSeparations',uniqueSeparations, ...
    'CombinedProducts',groupProducts, ...
    'DelayPerSeparation',groupDelays);
end
