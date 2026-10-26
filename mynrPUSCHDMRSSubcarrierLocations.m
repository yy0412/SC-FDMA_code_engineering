function [dmrssc,dshiftsnodata] = mynrPUSCHDMRSSubcarrierLocations(carrier,pusch)
%mynrPUSCHDMRSSubcarrierLocations Type 1/2 comb locations for custom PRBs.

nRBSC = double(carrier.SubcarriersPerPRB);
if pusch.DMRS.DMRSConfigurationType == 1
    % Type 1: comb-2 grid. For 12-SC PRBs this is [0 2 4 6 8 10].
    dmrssc = (0:2:nRBSC-1).';
    dshiftsnodata = 0:min(double(pusch.DMRS.NumCDMGroupsWithoutData),2)-1;
else
    % Type 2 keeps the standard two-adjacent-RE group pattern scaled by
    % the available PRB width.
    firstPair = [0 1].';
    secondPair = floor(nRBSC/2) + [0 1].';
    dmrssc = unique([firstPair; secondPair]);
    dmrssc = dmrssc(dmrssc >= 0 & dmrssc < nRBSC);
    dshiftsnodata = 2*(0:min(double(pusch.DMRS.NumCDMGroupsWithoutData),3)-1);
end
end
