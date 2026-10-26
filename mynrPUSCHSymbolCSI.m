function symbolCSI = mynrPUSCHSymbolCSI(carrier,pusch,reCSI)
%mynrPUSCHSymbolCSI Convert equalizer RE CSI to DFT-deprecoded symbol CSI.

msc = numel(pusch.PRBSet) * double(carrier.SubcarriersPerPRB);
reCSI = real(reCSI(:));

if pusch.TransformPrecoding
    symbolCSI = mynrTransformDeprecode(reCSI,msc) / sqrt(msc);
    symbolCSI = repmat(symbolCSI((1:msc:end).'),1,msc).';
    symbolCSI = symbolCSI(:);
else
    symbolCSI = reCSI;
end
symbolCSI = max(symbolCSI,eps);
end
