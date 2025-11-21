% +utils/reliability_split_half.m
function rel = reliability_split_half(Z_byOdor, num_splits, seed)
% split-half correlations across trials (time x cells frames flattened)
rel = nan(numel(Z_byOdor));
rng(seed);
for i = 1:numel(Z_byOdor)
    A = Z_byOdor{i};
    if isempty(A), continue; end
    for j = 1:numel(Z_byOdor)
        B = Z_byOdor{j};
        if isempty(B), continue; end
        cs = zeros(1,num_splits);
        for s = 1:num_splits
            idxA = randperm(size(A,3));
            idxB = randperm(size(B,3));
            A1   = mean(A(:,:,idxA(1:floor(end/2))),3);
            B1   = mean(B(:,:,idxB(1:floor(end/2))),3);
            cs(s)= corr(A1(:), B1(:), 'rows','complete');
        end
        rel(i,j) = mean(cs);
    end
end
end
