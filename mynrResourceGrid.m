function grid = mynrResourceGrid(carrier,numPorts)
%mynrResourceGrid Resource grid for a custom-width carrier.

if nargin < 2
    numPorts = 1;
end
grid = complex(zeros(carrier.NCustomSubcarriers,carrier.SymbolsPerSlot,numPorts));
end
