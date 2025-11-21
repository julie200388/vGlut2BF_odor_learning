function plot_traces_with_ci(Z_byOdor, resp, tvec, figdir, tag)
% RC = blue, NC = black, IC = green

colors = struct('RC',[0 0 1], 'NC',[0 0 0], 'IC',[0 0.5 0]);
fills  = struct('RC',[0.6 0.8 1], 'NC',[0.8 0.8 0.8], 'IC',[0.6 1 0.6]);

for k = 1:numel(Z_byOdor)
    Z = Z_byOdor{k};
    if isempty(Z), continue; end
    M = mean(Z, 3, 'omitnan');      % time x cells

    figure('Color','w'); hold on
    groups = {'RC','NC','IC'};              % fixed order
    hL = []; labels = {};                   % line handles + labels

    for gi = 1:numel(groups)
        gname = groups{gi};
        idx = resp(k).(gname);
        if isempty(idx), continue; end

        m  = mean(M(:,idx),2,'omitnan');
        se = std(M(:,idx),0,2,'omitnan')/sqrt(numel(idx));
        ci = 1.96*se;

        % shaded CI (hidden from legend)
        x = [tvec(:); flipud(tvec(:))];
        y = [m-ci; flipud(m+ci)];
        patch(x, y, fills.(gname), ...
              'EdgeColor','none', 'FaceAlpha',0.35, ...
              'HandleVisibility','off');

        % mean line (used in legend)
        h = plot(tvec, m, 'LineWidth',2, 'Color', colors.(gname), ...
                 'DisplayName', gname);
        hL(end+1) = h; %#ok<AGROW>
        labels{end+1} = sprintf('%s', gname); %#ok<AGROW>
    end

    xline(0,'k'); xline(2,'k');
    xlabel('Time (s)'); ylabel('zF/F');
    title(sprintf('Odor %d (%s)', k, tag));
    set(gca,'LineWidth',1.5,'FontSize',12); box off

    if ~isempty(hL)
        legend(hL, labels, 'Location','northeast', 'Box','off');
    end

    saveas(gcf, fullfile(figdir, sprintf('traces_%s_odor%02d.png',tag,k)));
    close
end
end
