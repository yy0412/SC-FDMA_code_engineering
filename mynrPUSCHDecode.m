function [cws,rxSymbols] = mynrPUSCHDecode(carrier,pusch,sym,nVar,varargin)
%mynrPUSCHDecode Decode custom-width PUSCH symbols to soft bits.
%
%   Based on nrPUSCHDecode. It uses the original toolbox functions for
%   MIMO deprecoding, demodulation, and descrambling. Only transform
%   deprecoding uses an explicit Msc.

if nargin < 4
    nVar = 1e-10;
end
validateConfig(pusch);

if isempty(pusch.NID)
    nid = carrier.NCellID;
else
    nid = pusch.NID(1);
end

if isempty(sym)
    deprecoded = zeros(0,pusch.NumLayers,'like',sym);
else
    if strcmpi(pusch.TransmissionScheme,'codebook')
        W = nrPUSCHCodebook(pusch.NumLayers,size(sym,2),pusch.TPMI,pusch.TransformPrecoding,pusch.CodebookType);
    else
        W = eye(pusch.NumLayers);
    end
    deprecoded = sym * pinv(W);
end

if pusch.TransformPrecoding
    msc = numel(pusch.PRBSet) * double(carrier.SubcarriersPerPRB);
    detransformed = mynrTransformDeprecode(deprecoded,msc);
else
    detransformed = deprecoded;
end

rxSymbols = detransformed;
demod = nrSymbolDemodulate(rxSymbols,pusch.Modulation,nVar,varargin{:});
cws = nrPUSCHDescramble(demod,nid,pusch.RNTI,pusch.NRAPID);
end
