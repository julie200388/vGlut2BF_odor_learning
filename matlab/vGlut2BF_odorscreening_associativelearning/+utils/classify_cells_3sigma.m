% +utils/classify_cells_3sigma.m
function resp = classify_cells_3sigma(MeanPerCell, FR, BS)
% outputs resp struct with fields RC/IC/NC (indices per odor)
resp = struct('RC',{[]}, 'IC',{[]}, 'NC',{[]});
resp = repmat(resp, 1, numel(MeanPerCell));
b0 = FR+1; b1 = BS*FR; r0 = BS*FR+1; r1 = (BS+2)*FR;
for k = 1:numel(MeanPerCell)
    M = MeanPerCell{k};
    if isempty(M), continue; end
    muB = mean(M(b0:b1,:), 1);
    sdB = std(M(b0:b1,:), 0, 1);
    sdB(sdB==0)=1;
    muR = mean(M(r0:r1,:), 1);
    RC = find(muR > 3*sdB);
    IC = find(muR < -3*sdB);
    allIdx = 1:size(M,2);
    NC = setdiff(allIdx, union(RC,IC));
    resp(k).RC = RC(:)'; resp(k).IC = IC(:)'; resp(k).NC = NC(:)';
end
end
