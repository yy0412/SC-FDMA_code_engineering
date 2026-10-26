function sym = mynrLowPAPRSequence(dmrs,nsc,prbset,nslot,nsymbol,ldash,symbperslot)
%mynrLowPAPRSequence Low-PAPR DM-RS sequence for custom PRB widths.
%
%   This is a narrow adaptation of nr5g.internal.pusch.lowPAPRSequence.
%   It keeps the same hopping calculation and calls nrLowPAPRS from the
%   5G Toolbox, which implements TS 38.211 Section 5.2.2.

if isempty(prbset)
    sym = complex([]);
    return;
end

[minprb,maxprb] = bounds(prbset);
mzc = nsc*(maxprb-minprb+1);
nsymbol = nsymbol - ldash;

if ~isfield(dmrs,'Type') || dmrs.Type == 1
    [u,v] = nr5g.internal.pusch.getHoppingParameters(dmrs,mzc,nslot,nsymbol,symbperslot);
    alpha = 0;
    base = localNrLowPAPRS(u(1),v(1),alpha,mzc);
else
    nidselected = min(dmrs.NSCID+1,length(dmrs.NIDNSCID));
    nidnscid = double(dmrs.NIDNSCID(nidselected));
    nscid = double(dmrs.NSCID);
    if isempty(dmrs.NRSID)
        dmrs.NRSID = nidnscid;
    end
    u = nr5g.internal.pusch.getHoppingParameters(dmrs,mzc,nslot,nsymbol,symbperslot);
    cinit = mod(2^17*(symbperslot*nslot + nsymbol + 1)*(2*nidnscid + 1) + 2*nidnscid + nscid,2^31);
    base = localNrLowPAPRS(u(1),cinit,mzc);
end

symmat = reshape(base,nsc,[]);
sym = reshape(symmat(:,prbset-minprb+1),[],1);
end

function seq = localNrLowPAPRS(a,b,c,m)
% nrLowPAPRS only defines short type-1 lengths 6/12/18/24 and otherwise
% ZC-derived lengths >= 30. For non-standard PRB widths that produce a
% shorter non-allowed length, generate the next valid TS sequence and take
% the scheduled part. Standard 12-SC PRBs are unaffected.
if nargin == 4
    mvalid = localValidLength(m);
    seq = nrLowPAPRS(a,b,c,mvalid);
else
    m = c;
    mvalid = localValidLength(m);
    seq = nrLowPAPRS(a,b,mvalid);
end
seq = seq(1:m);
end

function mvalid = localValidLength(m)
if m == 0 || m >= 30 || ismember(m,[6 12 18 24])
    mvalid = m;
elseif m < 6
    mvalid = 6;
elseif m < 12
    mvalid = 12;
elseif m < 18
    mvalid = 18;
elseif m < 24
    mvalid = 24;
else
    mvalid = 30;
end
end
