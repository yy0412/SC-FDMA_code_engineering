function grid = mynrOFDMDemodulate(carrier,waveform,varargin)
%mynrOFDMDemodulate OFDM demodulation without 12-subcarrier padding.
%
%   The FFT/IFFT scaling follows nrOFDMDemodulate so that the standard
%   N0 = 1/sqrt(Nfft*SNR) noise model maps to nVar = N0^2*Nfft per RE.

info = mynrOFDMInfo(carrier);
info.CyclicPrefixFraction = 0.5; % Used for cut the cyclic prefix 
[T,~] = size(waveform);
slotLength = sum(info.SymbolLengths);
if T < slotLength
    error('mynrOFDMDemodulate:WaveformTooShort','WAVEFORM is shorter than one custom OFDM slot.');
end

K = carrier.NCustomSubcarriers;
cpLengths = info.CyclicPrefixLengths;
firstSC = (info.Nfft/2) - (K/2) + 1;
nullIndices = [1:(firstSC-1) (firstSC+K):info.Nfft].';
symOffset = fix(cpLengths * info.CyclicPrefixFraction);

grid = ofdm_demodulate(waveform,info.Nfft,cpLengths,symOffset,nullIndices);

% Phase decompensation, TS 38.211 Section 5.4
if (any(info.SymbolPhases~=0))
    grid = grid .* exp(-1i * symbolPhases);
end


end

function y = ofdm_demodulate(x,nfft,cplen,symoffset,nullidx)

    p = struct();
    p.NumReceiveAntennas = size(x,2);
    p.FFTLength = nfft;
    p.NumSymbols = numel(cplen);
    p.NullIndices = nullidx;
    p.hasPilots = false;
    p.PilotIndices = [];
    p.SymbolOffset = symoffset;
    p.CyclicPrefixLength = cplen;
    y = comm.internal.ofdm.demodulate(x,p);

    dataIdx = true(nfft,1);
    dataIdx(nullidx) = false;
    y = y(dataIdx,:,:);

end