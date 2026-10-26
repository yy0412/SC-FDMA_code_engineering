function [waveform,info] = mynrOFDMModulate(carrier,grid,varargin)
%mynrOFDMModulate OFDM modulation without 12-subcarrier padding.
%
%   The FFT/IFFT scaling follows nrOFDMModulate so that the standard
%   N0 = 1/sqrt(Nfft*SNR) noise model maps to nVar = N0^2*Nfft per RE.

info = mynrOFDMInfo(carrier);
[K,L,P] = size(grid);
if K ~= carrier.NCustomSubcarriers || L ~= carrier.SymbolsPerSlot
    error('mynrOFDMModulate:GridSizeMismatch','GRID must be K-by-SymbolsPerSlot-by-P for the custom carrier.');
end

waveform = complex(zeros(sum(info.SymbolLengths),P,'like',grid));
active = info.ActiveSubcarrierIndices;

for p = 1:P % Leave an interface for antenna port
    offset = 0; % Leave an interface for delay
    for l = 1:L % For OFDM symbols
        f = complex(zeros(info.Nfft,1,'like',grid));
        f(active) = grid(:,l,p); % Allocate the active subcarreirs (near the central DC frequency)
        x = ifft(ifftshift(f),info.Nfft); % Move the negative frequecy band to the positive
        cp = x(end-info.CyclicPrefixLengths(l)+1:end);
        y = [cp; x]; 
        waveform(offset+(1:numel(y)),p) = y; % This is unwindowed waveform
        offset = offset + numel(y);
    end
    % Consider each OFDM symbol waveform edge problem, add the raised cosine window
    % Windowing (and optional overlapping)
    overlapping = true; 
    N = L;
    if ~isempty(info.Windowing)
        windowed = windowing(waveform,info,N,overlapping);
    else
        windowed = y;
    end
    waveform = windowed;
end
end

% This function is copied from nr5g.internal.OFDMModulate
function waveform = windowing(unwindowed,info,N,overlapping)

    % 'E' is the number of samples by which each OFDM symbol is cyclically
    % extended for the purposes of windowing and overlapping
    E = info.Windowing;
    
    % Calculate output waveform duration 'T' and number of antennas 'P'
    % and create empty waveform
    T = size(unwindowed,1);
    if (~overlapping)
        T = T + (E*N);
    end
    P = size(unwindowed,2);
    waveform = complex(zeros([T P],'like',unwindowed));
    
    % Initialize the cursors used to index the input and output waveforms
    % and the useful OFDM symbol period 'N_u'
    incursor = 0;
    outcursor = 0;
    N_u = info.Nfft;
    cpLengths = info.CyclicPrefixLengths;
    % For each OFDM symbol
    for n = 1:N
    
        % Get the cyclic prefix length 'N_CP' and total OFDM symbol
        % duration 'N_tot'
        N_CP = cpLengths(n);
        N_tot = N_CP + N_u;
    
        % Extract the current unwindowed OFDM symbol
        symbol = unwindowed(incursor + (1:N_tot),1:P);
        incursor = incursor + N_tot;
    
        % Create the window function (duration N_tot+E), cyclically extend
        % the OFDM symbol by a further 'E' samples at the start, and then
        % apply the window function (pointwise multiply)
        window = raised_cosine_window(N_tot,E);
        extended = [symbol(N_u-E + (1:E),1:P); symbol];
        windowed = extended .* window;
    
        % If overlapping is enabled and this is the last OFDM symbol
        if (overlapping && n==N)
    
            % Add the first 'N_tot' samples of the windowed OFDM symbol to
            % the output waveform
            outidx = outcursor + (1:N_tot);
            waveform(outidx,1:P) = waveform(outidx,1:P) + windowed(1:N_tot,1:P);
    
            % Add the last 'E' samples of the windowed OFDM symbol to
            % the first 'E' samples of the waveform
            waveform(1:E,1:P) = waveform(1:E,1:P) + windowed(end-E+1:end,1:P);
    
            % Perform a cyclic shift of -E samples to move the window
            % transition at the start of the waveform to the end of the
            % waveform. The first sample of the waveform is now the first
            % sample of the first cyclic prefix
            waveform = circshift(waveform,-E,1);
    
        else
    
            % Add the windowed OFDM symbol to the output waveform
            outidx = outcursor + (1:(N_tot + E));
            waveform(outidx,1:P) = waveform(outidx,1:P) + windowed;
    
        end
    
        % Update the output cursor to point to the start of where the next
        % OFDM symbol should be added to the waveform, taking into
        % consideration whether or not overlapping is enabled
        outcursor = outcursor + N_tot;
        if (~overlapping)
            outcursor = outcursor + E;
        end

    end

end

% Raised Cosine window creation; creates a window function of length n+e
% with raised cosine transitions on the first and last 'e' samples.
function p = raised_cosine_window(n,e)

    p = 0.5*(1-sin(pi*(e+1-2*(1:e).')/(2*e)));
    p = [p; ones(n-e,1); flipud(p)];

end

