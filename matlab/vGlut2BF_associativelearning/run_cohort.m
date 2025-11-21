function run_cohort(cohort)
% RUN_COHORT  Combine multiple animals for group figures.
% Usage:
cohort.animal_mat_files = { '/Users/peyshyuanchin/Documents/Matlabscriptforpaper/outputs/20251117_111645/data/103025_110725_H1_celltraces_animal_for_cohort.mat', '/Users/peyshyuanchin/Documents/Matlabscriptforpaper/outputs/20251117_115907/data/103025_110725_H2_celltraces_animal_for_cohort.mat'};
cohort.out_dir = './outputs_cohort';
cohort.FR = 15;   % frames/sec
cohort.BS = 3;    % baseline seconds before odor
%   run_cohort(cohort);
%%
assert(isfield(cohort,'animal_mat_files') && ~isempty(cohort.animal_mat_files), ...
    'Provide cohort.animal_mat_files (cellstr of .mat files).');

if ~isfield(cohort,'out_dir'), cohort.out_dir = './outputs_cohort'; end
dateTag = datestr(now,'yyyymmdd_HHMMSS');
outroot = fullfile(cohort.out_dir, dateTag);
figdir  = fullfile(outroot, 'figures'); 
datadir = fullfile(outroot, 'data');
utils.ensure_dir(outroot); utils.ensure_dir(figdir); utils.ensure_dir(datadir);

% ---------- LOAD ALL ANIMALS ----------
A = utils.load_animals_for_cohort(cohort.animal_mat_files);

% Sanity: assume all animals share same odor count within bf/af
K = min( cellfun(@numel,{A.MeanPerCell_bf}) );

% ---------- CONCATENATE CELLS ACROSS ANIMALS (for traces & heatmaps) ----------
% Build cohort-level mean-per-cell matrices by concatenating columns (cells)
C.MeanPerCell_bf = utils.concat_means_across_animals(A, K, 'bf');
C.MeanPerCell_af = utils.concat_means_across_animals(A, K, 'af');

K = min(numel(C.MeanPerCell_bf), numel(C.MeanPerCell_af));

% bf only, odors 2..K
utils.cohort_pca_and_plots({C.MeanPerCell_bf}, figdir, 'bf', ...
    'SetLabels', {'bf'}, 'OdorIdx', 2:K);

% af only, odors 2..K
utils.cohort_pca_and_plots({C.MeanPerCell_af}, figdir, 'af', ...
    'SetLabels', {'af'}, 'OdorIdx', 2:K);

% bf vs af overlaid, odors 2..K
utils.cohort_pca_and_plots({C.MeanPerCell_bf, C.MeanPerCell_af}, figdir, 'bfaf', ...
    'SetLabels', {'bf','af'}, 'OdorIdx', 2:K);

% Time axis (use the shortest across animals to be safe)
t_bf = utils.common_time(A, 'bf');
t_af = utils.common_time(A, 'af');

% ---------- PLOT: cohort average traces (across ALL cells) per odor ----------
utils.cohort_plot_allcells_traces(C.MeanPerCell_bf, C.MeanPerCell_af, t_bf, t_af, figdir);

% --- build reference order from odor 2 BEFORE condition ---
if numel(C.MeanPerCell_bf) >= 2 && ~isempty(C.MeanPerCell_bf{2})
    M2bf = C.MeanPerCell_bf{2};
    % mean response in [BS, BS+2] window is typical, but if you prefer a single ~2s index:
    resp2 = mean(M2bf(cohort.BS*cohort.FR+1:(cohort.BS+2)*cohort.FR, :), 1, 'omitnan');
    [~, order_ref] = sort(resp2, 'descend');
else
    warning('Odor 2 (bf) missing; falling back to odor 1 for sorting.');
    M1bf = C.MeanPerCell_bf{1};
    resp1 = mean(M1bf(cohort.BS*cohort.FR+1:(cohort.BS+2)*cohort.FR, :), 1, 'omitnan');
    [~, order_ref] = sort(resp1, 'descend');
end
% ---------- PLOT: cohort heatmaps (cells concatenated across animals) ----------
utils.cohort_plot_heatmaps( ...
    C.MeanPerCell_bf, ...
    C.MeanPerCell_af, ...
    t_bf, t_af, ...
    figdir, 'cohort');


% ---------- EUCLIDEAN DISTANCE (bf vs af): mean ± SD across animals ----------
D = utils.cohort_euclidean_distance(A, K);  % struct with fields: t, mean, sd, per-odor
save(fullfile(datadir,'euclidean_distance_bf_vs_af.mat'), 'D');

% Plot mean ± SD for each odor
for k = 1:K
    figure('Color','w'); hold on
    plot(D(k).t, D(k).mean, 'LineWidth', 2, 'Color', [0 0 0]);
    % SD band (not CI): mean ± SD across animals
    x = [D(k).t(:); flipud(D(k).t(:))];
    y = [D(k).mean - D(k).sd; flipud(D(k).mean + D(k).sd)];
    patch(x,y,[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha', 0.5);
    uistack(findobj(gca,'Type','line','-and','Color',[0 0 0]),'top');
    xline(0,'k'); xline(2,'k');
    xlabel('Time (s)'); ylabel('Euclidean distance (bf vs af)');
    title(sprintf('Cohort Euclidean distance — Odor %d (mean ± SD across animals)', k));
    set(gca,'LineWidth',1.5,'FontSize',12); box off
    saveas(gcf, fullfile(figdir, sprintf('cohort_euclid_odor%02d.png',k)));
    close
end

fprintf('✓ Cohort plots saved to %s\n', figdir);
%%---------- Crosstab (AF): Odor A vs Odor B organized by A-status ----------
FR = cohort.FR;  % or set explicitly, e.g., 15
BS = cohort.BS;  % e.g., 3
odors_pair = [2 3];  % Pentanol (2) vs Hexanol (3)

XT = utils.cohort_af_crosstab(A, FR, BS, odors_pair);

% Save cohort summary tables
cohort_counts   = array2table(XT.cohort.counts,  'VariableNames', XT.labels.cols, 'RowNames', XT.labels.rows);
cohort_rowperc  = array2table(XT.cohort.rowperc, 'VariableNames', XT.labels.cols, 'RowNames', XT.labels.rows);
writetable(cohort_counts,  fullfile(datadir, 'crosstab_counts_AF_odorA_vs_odorB.csv'), 'WriteRowNames', true);
writetable(cohort_rowperc, fullfile(datadir, 'crosstab_rowperc_AF_odorA_vs_odorB.csv'), 'WriteRowNames', true);

% Optional: per-animal export (one tidy CSV)
rows = {};
for a = 1:numel(XT.per_animal)
    C = XT.per_animal(a).counts;
    RP = XT.per_animal(a).rowperc;
    for r = 1:3
        rows(end+1,:) = {string(XT.per_animal(a).animal_id), XT.labels.rows{r}, 'B_RC', C(r,1), RP(r,1)}; %#ok<AGROW>
        rows(end+1,:) = {string(XT.per_animal(a).animal_id), XT.labels.rows{r}, 'B_IC', C(r,2), RP(r,2)}; %#ok<AGROW>
        rows(end+1,:) = {string(XT.per_animal(a).animal_id), XT.labels.rows{r}, 'B_NC', C(r,3), RP(r,3)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'animal_id','A_status','B_status','count','row_percent'});
writetable(T, fullfile(datadir, 'crosstab_per_animal_AF_odorA_vs_odorB.csv'));

end
