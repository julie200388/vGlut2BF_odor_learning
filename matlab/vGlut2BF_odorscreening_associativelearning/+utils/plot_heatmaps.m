% +utils/plot_heatmaps.m
function plot_heatmaps(Z_byOdor, tvec, FR, BS, figdir, tag)
TF = numel(tvec);
for k = 1:numel(Z_byOdor)
    Z = Z_byOdor{k};
    if isempty(Z), continue; end
    M = mean(Z,3,'omitnan');  % time x cells
    % optional: sort by response at (BS+2)s mark
    [~,ord] = sort(M((BS+2)*FR, :), 'descend');
    M = M(:,ord);
    figure('Color','w'); 
    h = heatmap(M','Colormap',parula,'GridVisible','off'); %#ok<NASGU>
    caxis([-1 2]);
    % vertical lines at odor on/off
    S = struct(h); ax = S.Axes;
    xline(ax, find(tvec==0)+0.5, 'k', 'LineWidth',2);
    xline(ax, find(tvec==2)+0.5, 'k', 'LineWidth',2);
    title(sprintf('Heatmap Odor %d (%s)',k,tag));
    saveas(gcf, fullfile(figdir, sprintf('heatmap_%s_odor%02d.png',tag,k)));
    close
end
end
