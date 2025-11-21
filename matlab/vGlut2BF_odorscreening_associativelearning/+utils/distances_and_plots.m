% +utils/distances_and_plots.m
function distances_and_plots(Z_byOdor_bf, Z_byOdor_af, tvec_bf, tvec_af, figdir)
% Example: Euclidean distance over time between bf vs af per odor
K = min(numel(Z_byOdor_bf), numel(Z_byOdor_af));
for k = 1:K
    A = Z_byOdor_bf{k}; B = Z_byOdor_af{k};
    if isempty(A) || isempty(B), continue; end
    mA = mean(A,3,'omitnan'); % time x cells
    mB = mean(B,3,'omitnan');
    T  = min(size(mA,1), size(mB,1));
    d  = sqrt(sum((mA(1:T,:) - mB(1:T,:)).^2, 2));
    figure('Color','w'); 
    plot( tvec_bf(1:T), d, 'LineWidth',2 );
    xlabel('Time (s)'); ylabel('Euclidean distance (bf vs af)');
    title(sprintf('Odor %d distance',k));
    xline(0); xline(2);
    box off; set(gca,'LineWidth',1.5,'FontSize',12);
    saveas(gcf, fullfile(figdir, sprintf('dist_bf_vs_af_odor%02d.png',k)));
    close
end
end
