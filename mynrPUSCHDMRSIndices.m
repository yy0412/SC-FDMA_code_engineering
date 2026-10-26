function ind = mynrPUSCHDMRSIndices(carrier,pusch,varargin)
%mynrPUSCHDMRSIndices PUSCH DM-RS RE indices for custom-width PRBs.
%
%   Based on nrPUSCHDMRSIndices. Type 1 uses a comb-2 grid; for a normal
%   12-subcarrier PRB this exactly reduces to the toolbox pattern
%   [0 2 4 6 8 10] plus the port delta shift.

validateConfig(pusch);
nRBSC = double(carrier.SubcarriersPerPRB);
nSizeGrid = double(carrier.NSizeGrid);
symbperslot = double(carrier.SymbolsPerSlot);
nTotalRE = nRBSC * nSizeGrid * symbperslot;

[~,dmrssymbols] = mynrPUSCHDMRSSymbols(carrier,pusch);
prbset = reshape(double(pusch.PRBSet),1,[]);
[dmrssc,~] = mynrPUSCHDMRSSubcarrierLocations(carrier,pusch);

if isempty(pusch.DMRS.DMRSPortSet)
    layers = 0:double(pusch.NumLayers)-1;
else
    layers = double(pusch.DMRS.DMRSPortSet);
end
nLayers = numel(layers);
if strcmpi(pusch.TransmissionScheme,'codebook')
    nports = double(pusch.NumAntennaPorts);
else
    nports = nLayers;
end

% Delta shifts follow the standard config object where available, but are
% wrapped into the custom PRB width.
dshifts = double(pusch.DMRS.DeltaShifts(:,1:nLayers));
dshifts = mod(dshifts,nRBSC);

indPerLayer = zeros(0,nLayers);
for i = 1:numel(dmrssymbols)
    symI = dmrssymbols(i);
    for prb = prbset
        base = nRBSC*prb + nRBSC*nSizeGrid*symI;
        k = repmat(dmrssc,1,nLayers) + repmat(dshifts,size(dmrssc,1),1);
        k = mod(k,nRBSC);
        indPerLayer = [indPerLayer; k + base + repmat(nTotalRE*(0:nLayers-1),size(k,1),1)]; %#ok<AGROW>
    end
end

if strcmpi(pusch.TransmissionScheme,'codebook') && ~isempty(indPerLayer)
    planeoffsets = nTotalRE*(0:nports-1);
    ind = repmat(indPerLayer(:,1),1,nports) + repmat(planeoffsets,size(indPerLayer,1),1);
else
    ind = indPerLayer;
end

if nargin > 2
    opts = nr5g.internal.parseOptions('mynrPUSCHDMRSIndices',{'IndexStyle','IndexBase'},varargin{:});
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
