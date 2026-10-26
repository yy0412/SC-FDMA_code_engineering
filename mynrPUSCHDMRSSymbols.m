function [symbolset,dmrssymbols,ldash] = mynrPUSCHDMRSSymbols(carrier,pusch)
%mynrPUSCHDMRSSymbols PUSCH and DM-RS OFDM symbol sets.
%
%   This follows the same lookup used by nrPUSCHDMRS/nrPUSCHIndices via
%   nr5g.internal.pusch.lookupPUSCHDMRSSymbols, with a compact fallback for
%   custom-symbol-set operation.

l0 = pusch.SymbolAllocation(1);
lalloc = pusch.SymbolAllocation(2);
symbolset = l0 + (0:lalloc-1);

if ~isempty(pusch.DMRS.CustomSymbolSet)
    dmrssymbols = reshape(double(pusch.DMRS.CustomSymbolSet),1,[]);
    dmrssymbols = dmrssymbols(ismember(dmrssymbols,symbolset));
else
    ftable.ChannelName = 'PUSCH';
    ftable.MappingTypeB = strcmpi(pusch.MappingType,'B');
    ftable.DMRSSymbolSet = @nr5g.internal.pusch.lookupPUSCHDMRSSymbols;
    ftable.IntraSlotFreqHoppingFlag = strcmpi(pusch.FrequencyHopping,'intraSlot');
    nSizeBWP = pusch.NSizeBWP;
    if isempty(nSizeBWP)
        nSizeBWP = carrier.NSizeGrid;
    end
    [~,symbolset,dmrssymbols,ldash] = nr5g.internal.pxsch.initializeResources(carrier,pusch,nSizeBWP,ftable);
    return;
end

ldash = zeros(1,numel(dmrssymbols));
if pusch.DMRS.DMRSLength == 2
    for i = 2:numel(dmrssymbols)
        if any(dmrssymbols(i) == dmrssymbols(1:i-1)+1)
            ldash(i) = 1;
        end
    end
end
end
