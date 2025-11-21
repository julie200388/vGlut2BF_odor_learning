function cohort_plot_allcells_traces(MeanPerCell_bf, MeanPerCell_af, t_bf, t_af, figdir)
% Average across ALL concatenated cells, per odor (bf=black, af=red), with 95% CI across cells
K = min(numel(MeanPerCell_bf), numel(MeanPerCell_af));
for k = 1:K
    Mb = MeanPerCell_bf{k}; 
    Ma = MeanPerCell_af{k}; 
    if isempty(Mb) || isempty(Ma), continue; end
    T = min(size(Mb,1), size(Ma,1));
    t = t_bf(1:T);

    [mb, lb, hb] = mean_ci95_cells(Mb(1:T,:));
    [ma, la, ha] = mean_ci95_cells(Ma(1:T,:));

    figure('Color','w'); hold on
    plot(t, mb, 'k','LineWidth',2);
    patch([t;flipud(t)],[lb;flipud(hb)], [0.8 0.8 0.8],'EdgeColor','none','FaceAlpha', 0.5);
    plot(t, ma, 'Color',[0.8 0 0],'LineWidth',2);
    patch([t;flipud(t)],[la;flipud(ha)], [0.9 0.6 0.6],'EdgeColor','none','FaceAlpha', 0.5);
    uistack(findobj(gca,'Type','line','-and','Color',[0 0 0]),'top');
    uistack(findobj(gca,'Type','line','-and','Color',[0.8 0 0]),'top');
    xline(0,'k'); xline(2,'k');
    xlabel('Time (s)'); ylabel('zF/F');
    title(sprintf('Cohort average across ALL cells — Odor %d', k));
    set(gca,'LineWidth',1.5,'FontSize',12); box off
    saveas(gcf, fullfile(figdir, sprintf('cohort_allcells_overlay_odor%02d.png',k)));
    close
end
end

function [m, lo, hi] = mean_ci95_cells(X)
m  = mean(X,2,'omitnan');
se = std(X,0,2,'omitnan') ./ max(1, sqrt(sum(~isnan(X),2)));
ci = 1.96 * se;
lo = m - ci; hi = m + ci;
end
