function plot_allcells_overlay(MeanPerCell_A, MeanPerCell_B, t_A, t_B, figdir, nameTag)
% PLOT_ALLCELLS_OVERLAY
% Overlay mean + 95% CI across all neurons for condition A vs condition B.
% nameTag = string used in filenames so outputs don't overwrite.

K = min(numel(MeanPerCell_A), numel(MeanPerCell_B));

for k = 1:K
    MA = MeanPerCell_A{k};
    MB = MeanPerCell_B{k};
    if isempty(MA) || isempty(MB), continue; end

    % Match time lengths
    T = min(size(MA,1), size(MB,1));
    t  = t_A(1:T);
    MA = MA(1:T,:);
    MB = MB(1:T,:);

    % Compute mean + CI
    [mA, loA, hiA] = mean_ci_95(MA);
    [mB, loB, hiB] = mean_ci_95(MB);

    % ----- Plot -----
    figure('Color','w'); hold on; t = t(:);

    % A: black
    plot(t, mA, 'k', 'LineWidth', 2);
    patch([t; flipud(t)], [loA; flipud(hiA)], [0.8 0.8 0.8], ...
          'EdgeColor','none','FaceAlpha',0.5);

    % B: red
    plot(t, mB, 'Color',[0.8 0 0], 'LineWidth',2);
    patch([t; flipud(t)], [loB; flipud(hiB)], [0.9 0.6 0.6], ...
          'EdgeColor','none','FaceAlpha',0.5);

    % Make sure lines are on top
    uistack(findobj(gca,'Type','line','-and','Color',[0 0 0]),'top');
    uistack(findobj(gca,'Type','line','-and','Color',[0.8 0 0]),'top');

    xline(0,'k'); xline(2,'k');
    xlabel('Time (s)'); ylabel('zF/F');
    title(sprintf('%s — Odor %d', strrep(nameTag,'_',' '), k));
    set(gca,'LineWidth',1.5,'FontSize',12); box off;

    % ---------- Unique filename ----------
    fname = sprintf('allcells_overlay_%s_odor%02d.png', nameTag, k);
    saveas(gcf, fullfile(figdir, fname));
    close
end
end


function [m, lo, hi] = mean_ci_95(X)
% X: [time x cells]
m  = mean(X,2,'omitnan');
se = std(X,0,2,'omitnan') ./ max(1, sqrt(sum(~isnan(X),2)));
ci = 1.96 * se;
lo = m - ci;
hi = m + ci;
end
