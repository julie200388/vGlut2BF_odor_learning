function cohort_plot_heatmaps(MeanPerCell_bf, MeanPerCell_af, t_bf, t_af, figdir, tag)
% COHORT_PLOT_HEATMAPS
% One figure per odor:
%   left  = BEFORE heatmap
%   right = AFTER heatmap
% Cells are sorted based on odor 2 (before condition).

K = numel(MeanPerCell_bf);

% --- find index near 0 s and 2 s for bf/af ---
[~, i0_bf] = min(abs(t_bf - 0));
[~, i2_bf] = min(abs(t_bf - 2));
[~, i0_af] = min(abs(t_af - 0));
[~, i2_af] = min(abs(t_af - 2));

% --- build reference order from odor 2 BEFORE condition ---
if K >= 2 && ~isempty(MeanPerCell_bf{2})
    Mref = MeanPerCell_bf{2};
    % sort by value near 2 s
    [~, idx2] = min(abs(t_bf - 2));
    [~, order_ref] = sort(Mref(idx2, :), 'descend');
else
    warning('Odor 2 (bf) missing or empty; using odor 1 as reference for sorting.');
    Mref = MeanPerCell_bf{1};
    [~, idx2] = min(abs(t_bf - 2));
    [~, order_ref] = sort(Mref(idx2, :), 'descend');
end

for k = 1:K
    Mb = MeanPerCell_bf{k};
    Ma = MeanPerCell_af{k};
    if isempty(Mb) || isempty(Ma), continue; end

    % apply same cell order (if dimensions allow)
    if max(order_ref) <= size(Mb,2)
        Hb = Mb(:, order_ref);
    else
        Hb = Mb;
        warning('Odor %d (bf) has fewer cells than reference sorting order.', k);
    end

    if max(order_ref) <= size(Ma,2)
        Ha = Ma(:, order_ref);
    else
        Ha = Ma;
        warning('Odor %d (af) has fewer cells than reference sorting order.', k);
    end

    % -------- FIGURE --------
    fig = figure('Color','w', 'Position', [300 300 1200 500]);

    % BEFORE subplot
    subplot(1,2,1);
    h1 = heatmap(Hb', 'Colormap', parula, 'GridVisible', 'off');
    caxis([-1 2]);
    S1 = struct(h1); ax1 = S1.Axes;
    xline(ax1, i0_bf+0.5, 'k', 'LineWidth', 2);
    xline(ax1, i2_bf+0.5, 'k', 'LineWidth', 2);
    title(sprintf('Before — Odor %d', k));

    % AFTER subplot
    subplot(1,2,2);
    h2 = heatmap(Ha', 'Colormap', parula, 'GridVisible', 'off');
    caxis([-1 2]);
    S2 = struct(h2); ax2 = S2.Axes;
    xline(ax2, i0_af+0.5, 'k', 'LineWidth', 2);
    xline(ax2, i2_af+0.5, 'k', 'LineWidth', 2);
    title(sprintf('After — Odor %d', k));

    % shared title + save
    sgtitle(sprintf('Cohort heatmap — Odor %d (%s, sorted by odor 2 bf)', k, tag), ...
            'FontSize', 16, 'FontWeight', 'bold');

    saveas(fig, fullfile(figdir, sprintf('cohort_heatmap_pair_%s_odor%02d.png', tag, k)));
    close(fig);
end
end


