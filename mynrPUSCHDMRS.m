function sym = mynrPUSCHDMRS(carrier,pusch,varargin)
%mynrPUSCHDMRS PUSCH DM-RS symbols for custom-width PRBs.
%
%   Based on nrPUSCHDMRS. For transform precoding, this function keeps
%   the 5G Toolbox low-PAPR sequence chain:
%   getHoppingParameters -> nrLowPAPRS (TS 38.211 Section 5.2.2).

validateConfig(pusch);
[~,dmrssymbols,ldash] = mynrPUSCHDMRSSymbols(carrier,pusch);
[dmrssc,~] = mynrPUSCHDMRSSubcarrierLocations(carrier,pusch);
ndmrsre = numel(dmrssc);
prbset = reshape(double(pusch.PRBSet),1,[]);
symbperslot = double(carrier.SymbolsPerSlot);
nslot = mod(double(carrier.NSlot),carrier.SlotsPerFrame);
nLayers = double(pusch.NumLayers);

if isempty(dmrssymbols)
    dmrsSymPort = complex(zeros(0,nLayers));
else
    if pusch.TransformPrecoding
        if iscell(pusch.Modulation)
            modstring = pusch.Modulation{1};
        else
            modstring = pusch.Modulation;
        end
        lowpaprtype = 1 + (pusch.DMRS.DMRSUplinkTransformPrecodingR16 && strcmpi(modstring,'pi/2-BPSK'));
        if isempty(pusch.DMRS.NRSID) && lowpaprtype == 1
            nrsid = carrier.NCellID;
        else
            nrsid = pusch.DMRS.NRSID;
        end
        nscid = pusch.DMRS.NSCID;
        if isempty(pusch.DMRS.NIDNSCID)
            nidnscid = carrier.NCellID;
        else
            nidnscid = pusch.DMRS.NIDNSCID;
        end
        dmrsStruct = struct('Type',lowpaprtype, ...
            'NIDNSCID',nidnscid,'NSCID',nscid(1), ...
            'NRSID',nrsid, ...
            'GroupHopping',pusch.DMRS.GroupHopping, ...
            'SequenceHopping',pusch.DMRS.SequenceHopping);
    else
        error('mynrPUSCHDMRS:CPODMUnsupported', ...
            'This custom PRB-width DM-RS implementation is intended for TransformPrecoding=true.');
    end

    cnLen = zeros(1,numel(dmrssymbols));
    symcell = cell(1,numel(dmrssymbols));
    for i = 1:numel(dmrssymbols)
        symcell{i} = mynrLowPAPRSequence(dmrsStruct,ndmrsre,prbset,nslot,dmrssymbols(i),ldash(i),symbperslot);
        cnLen(i) = numel(symcell{i});
    end

    cndmrs = [0 cumsum(cnLen)];
    dmrsSymPort = complex(zeros(cndmrs(end),nLayers));
    fmaskAllPorts = pusch.DMRS.FrequencyWeights;
    tmaskAllPorts = pusch.DMRS.TimeWeights;
    nwf = size(fmaskAllPorts,1);
    nfmaskrep = lcm(max(ndmrsre,1),nwf)/nwf;

    for pidx = 1:nLayers
        fmask = fmaskAllPorts(:,pidx);
        fmaskprb = reshape(repmat(fmask,nfmaskrep,1),ndmrsre,[]);
        tmask = tmaskAllPorts(ldash+1,pidx);
        for i = 1:numel(dmrssymbols)
            rbidxs = prbset;
            fmaskallprb = reshape(fmaskprb(:,mod(rbidxs,size(fmaskprb,2))+1),[],1);
            dmrsSymPort(cndmrs(i)+1:cndmrs(i+1),pidx) = symcell{i} .* fmaskallprb * tmask(i);
        end
    end
end

if strcmpi(pusch.TransmissionScheme,'codebook') && ~isempty(dmrsSymPort)
    W = nrPUSCHCodebook(nLayers,pusch.NumAntennaPorts,pusch.TPMI,pusch.TransformPrecoding,pusch.CodebookType);
    dmrs = dmrsSymPort * W;
else
    dmrs = dmrsSymPort;
end

if nargin > 2
    opts = nr5g.internal.parseOptions('mynrPUSCHDMRS',{'OutputDataType'},varargin{:});
    sym = cast(dmrs,opts.OutputDataType);
else
    sym = double(dmrs);
end
end
