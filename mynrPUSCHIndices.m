function [ind,info] = mynrPUSCHIndices(carrier,pusch,varargin)
%mynrPUSCHIndices PUSCH data RE indices with custom subcarriers per PRB.
%
%   Based on nrPUSCHIndices. The key change is replacing nRBSC = 12 by
%   carrier.SubcarriersPerPRB.

validateConfig(pusch);
nRBSC = double(carrier.SubcarriersPerPRB);    % Number of subcarriers per PRB
nSizeGrid = double(carrier.NSizeGrid);        % Number of PRBs
symbperslot = double(carrier.SymbolsPerSlot); % Number of OFDM symbols per slot
nTotalRE = nRBSC * nSizeGrid * symbperslot;   % Number of total REs per slot

[symbolset,dmrssymbols] = mynrPUSCHDMRSSymbols(carrier,pusch); % Generate the 0-based indices set for  
                                                               % allocated OFDM symbols and dmrs
prbset = reshape(double(pusch.PRBSet),1,[]); % Generate the 0-based indices set for PRB

% Define a new function to generate the 0-based dmrs indices for one basic PRB and 
% Delta shift for CDM groups without data (step size set in frequency domain) 
[dmrssc,dshiftsnodata] = mynrPUSCHDMRSSubcarrierLocations(carrier,pusch); % Depend on the custom subcarriers per PRB

% Non DM-RS resource elements in a DM-RS containing symbol
% Due to enabled transform precoding, the PUSCH data is banned within a DM-RS symbol
fullprb = true(nRBSC,1);
dmrsre = dmrssc + reshape(dshiftsnodata,1,[]);
dmrsre = unique(dmrsre(:));
dmrsre = dmrsre(dmrsre >= 0 & dmrsre < nRBSC);
fullprb(dmrsre+1) = false;
puschre = find(fullprb)-1; % Find PUSCH (non DM-RS) RE in a DM-RS containing symbol

% Calculate the PUSCH data RE indices
slotindices = zeros(0,1);
for i = 1:numel(symbolset) % Loop for each OFDM symbol
    symI = symbolset(i);
    if ismember(symI,dmrssymbols) % Decide if it is DM-RS symbol location
        recell = puschre;
    else
        recell = (0:nRBSC-1).';   % Allocate 0-based PUSCH RE indices to the current symbol
    end
    for prb = prbset % Loop for each PRB
        k0 = nRBSC*prb; % First subcarrier index for the current PRB
        slotindices = [slotindices; recell + k0 + nRBSC*nSizeGrid*symI]; %#ok<AGROW>
    end
end

% Capture number of ports (Currently only one antenna port is used)
if strcmpi(pusch.TransmissionScheme,'codebook')
    nports = double(pusch.NumAntennaPorts);
else
    nports = double(pusch.NumLayers);
end
ind = repmat(slotindices,1,nports) + repmat(nTotalRE*(0:nports-1),numel(slotindices),1);

% Modulation order for each codeword and combine information into a structure
qm = nr5g.internal.getQm(pusch.Modulation);
info.G = numel(slotindices) * qm * double(pusch.NumLayers); % Bit capacity of the PUSCH. This must be the
                                                            % length of codeword from the UL-SCH transport channel
info.Gd = numel(slotindices); % Number of REs per layer/port (symbol capacity of the PUSCH)
info.NREPerPRB = nRBSC*(numel(symbolset)-numel(dmrssymbols)) + numel(puschre)*numel(dmrssymbols); % Number of RE per PRB allocated to PUSCH (not
                                                                                                  % accounting for any reserved resources)
info.DMRSSymbolSet = dmrssymbols; % The OFDM symbol locations in a slot containing DM-RS (0-based)
info.PTRSSymbolSet = zeros(1,0);
info.PRBSet = prbset; % PRBs allocated for PUSCH within the BWP
info.SubcarriersPerPRB = nRBSC;

% Apply options
if nargin > 2
    opts = nr5g.internal.parseOptions('mynrPUSCHIndices',{'IndexStyle','IndexBase'},varargin{:});
    ind = localApplyIndexOptions(ind,[nRBSC*nSizeGrid symbperslot nports],opts);
else
    ind = uint32(ind + 1);
end
end

function out = localApplyIndexOptions(ind,dims,opts)
if strcmpi(opts.IndexBase,'1based')
    ind1 = ind + 1;
else
    ind1 = ind;
end
if strcmpi(opts.IndexStyle,'subscript')
    [k,l,p] = ind2sub(dims,ind+1);
    if strcmpi(opts.IndexBase,'0based')
        k = k-1; l = l-1; p = p-1;
    end
    out = [k(:) l(:) p(:)];
else
    out = ind1;
end
end
