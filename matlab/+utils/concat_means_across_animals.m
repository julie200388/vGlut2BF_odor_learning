function MeanPerCell_cohort = concat_means_across_animals(A, K, which)
% Concatenate time×cells mean-per-cell across animals (columns)
% which: 'bf' or 'af'
MeanPerCell_cohort = cell(1,K);
for k = 1:K
    Mlist = {};
    for a = 1:numel(A)
        if strcmp(which,'bf')
            M = A(a).MeanPerCell_bf{k};
            t = A(a).tvec_bf;
        else
            M = A(a).MeanPerCell_af{k};
            t = A(a).tvec_af;
        end
        if isempty(M), continue; end
        % Align to minimal common time length across animals for this odor
        if isempty(Mlist)
            T = size(M,1);
            tref = t(:);
        else
            T = min(T, size(M,1));
        end
        Mlist{end+1} = M; %#ok<AGROW>
    end
    % trim to T and concat columns
    X = [];
    for i = 1:numel(Mlist)
        X = [X, Mlist{i}(1:T,:)]; %#ok<AGROW>
    end
    MeanPerCell_cohort{k} = X; % time x cells_all
end
end
