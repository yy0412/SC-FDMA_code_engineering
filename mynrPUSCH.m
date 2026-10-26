function sym = mynrPUSCH(carrier,pusch,cw,varargin)
%mynrPUSCH PUSCH symbols for custom-width PRBs.
%
%   This is based on nrPUSCH and uses the original toolbox functions for
%   scrambling, modulation, layer mapping, and codebook precoding. Only
%   transform precoding uses mynrTransformPrecode because nrTransformPrecode
%   hardcodes Msc = MRB*12.

validateConfig(pusch);
if isempty(pusch.NID)
    nid = carrier.NCellID;
else
    nid = pusch.NID(1);
end

mods = nr5g.internal.validatePXSCHModulation('mynrPUSCH',pusch.Modulation,pusch.NumCodewords, ...
    {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'});
scrambled = nrPUSCHScramble(cw,nid,pusch.RNTI,pusch.NRAPID);
modulated = nrSymbolModulate(scrambled,mods{1},varargin{:});
layered = nrLayerMap(modulated,pusch.NumLayers);

if pusch.TransformPrecoding
    msc = numel(pusch.PRBSet) * double(carrier.SubcarriersPerPRB);
    transformed = mynrTransformPrecode(layered,msc);
else
    transformed = layered;
end

if strcmpi(pusch.TransmissionScheme,'codebook')
    W = nrPUSCHCodebook(pusch.NumLayers,pusch.NumAntennaPorts,pusch.TPMI,pusch.TransformPrecoding,pusch.CodebookType);
    sym = transformed * W;
else
    sym = transformed;
end
end
