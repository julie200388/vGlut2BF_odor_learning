function D = cohort_euclidean_distance(A, K)
% For each animal and odor:
%   compute Euclidean distance over time between mean bf vs mean af (time x cells)
% Then aggregate across animals: D(k).t, D(k).mean, D(k).sd

D = struct('t',{},'mean',{},'sd',{});
for k = 1:K
    % collect per-animal distances and time vectors
    dist_mat = [];  % time x animals
    t_common = [];
    for a = 1:numel(A)
        Mb = A(a).MeanPerCell_bf{k};
        Ma = A(a).MeanPerCell_af{k};
        if isempty(Mb) || isempty(Ma), continue; end
        T = min(size(Mb,1), size(Ma,1));
        t  = A(a).tvec_bf(1:T);

        db = sqrt(sum((Mb(1:T,:) - Ma(1:T,:)).^2, 2)); % Euclid per time
        dist_mat = [dist_mat, db]; %#ok<AGROW>

        if isempty(t_common)
            t_common = t(:);
        else
            % align to minimal length
            L = min(numel(t_common), numel(t));
            t_common = t_common(1:L);
            dist_mat = dist_mat(1:L, :);
        end
    end
    D(k).t    = t_common;
    D(k).mean = mean(dist_mat, 2, 'omitnan');
    D(k).sd   = std(dist_mat, 0, 2, 'omitnan');   % SD across animals
end
end
