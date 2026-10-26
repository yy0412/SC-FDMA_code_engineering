function info = mynrOFDMInfo(carrier)
%mynrOFDMInfo OFDM information for the custom subcarrier grid.
%
%   The time-domain OFDM parameters follow nrOFDMInfo for the configured
%   numerology. Only the number and placement of active subcarriers is
%   changed to K = NSizeGrid*SubcarriersPerPRB.

K = double(carrier.NCustomSubcarriers);

stdInfo = nrOFDMInfo(carrier);
if K > stdInfo.Nfft
    stdInfo = nrOFDMInfo(carrier,'Nfft',2^nextpow2(K));
end

symbolsPerSlot = double(stdInfo.SymbolsPerSlot);
slotInSubframe = mod(double(carrier.NSlot),double(stdInfo.SlotsPerSubframe));
symbolStart = slotInSubframe*symbolsPerSlot + 1;
symbolSet = symbolStart:(symbolStart + symbolsPerSlot - 1);

cpAll = stdInfo.CyclicPrefixLengths(:);
if numel(cpAll) >= symbolSet(end)
    cpLengths = cpAll(symbolSet);
else
    cpLengths = cpAll(1:symbolsPerSlot);
end

if isfield(stdInfo,'SymbolPhases')
    phaseAll = stdInfo.SymbolPhases(:);
    if numel(phaseAll) >= symbolSet(end)
        symbolPhases = phaseAll(symbolSet);
    else
        symbolPhases = phaseAll(1:symbolsPerSlot);
    end
else
    symbolPhases = zeros(symbolsPerSlot,1);
end

firstActive = floor((stdInfo.Nfft-K)/2) + 1;
active = firstActive:firstActive+K-1; % This is used for OFDM modulation (find the central (DC) frequency location)

info = stdInfo;
info.CyclicPrefixLengths = cpLengths;
info.SymbolLengths = cpLengths + stdInfo.Nfft;
info.SymbolPhases = symbolPhases;
info.NCustomSubcarriers = K;
info.ActiveSubcarrierIndices = active(:);
end
