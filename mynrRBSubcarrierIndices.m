function k = mynrRBSubcarrierIndices(carrier,prbset)
%mynrRBSubcarrierIndices 1-based subcarrier indices for custom-width PRBs.

nRBSC = double(carrier.SubcarriersPerPRB);
prbset = reshape(double(prbset),1,[]);
k = zeros(nRBSC*numel(prbset),1);
pos = 1;
for prb = prbset
    k(pos:pos+nRBSC-1) = prb*nRBSC + (1:nRBSC);
    pos = pos + nRBSC;
end
end
