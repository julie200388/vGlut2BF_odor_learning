% +utils/plot_heatmaps_multiday.m
function plot_heatmaps_multiday(Z_bf, t_bf, Z_af, t_af, FR, BS, figdir, tag, Z_day3, t_day3)
% PLOT_HEATMAPS_MULTIDAY
%   - Sorts cells according to response to **odor 2, day 1 (bf)**.
%   - For each odor k, plots bf / af / (optional) day3 heatmaps
%     in subplots of the same figure.
%
% Inputs:
%   Z_bf, Z_af, Z_day3 : 1xK cell arrays, each {k} = time x cells x trials
%   t_bf, t_af, t_day3 : time vectors for each day
%   FR, BS             : frame rate and baseline seconds
%   figdir             : output directory
%   tag                : string for filenames (e.g. 'single' or 'cohort')
%   Z_day3, t_day3     : optional (can be {} / [] if no third day)

if nargin < 9
    Z_day3 = {};
    t_day3 = [];
end

% ---------- determine how many "days" we actually have ----------
days = struct([]);
days(1).Z     = Z_bf;
days(1).t     = t_bf;
days(1).label = 'Day 1 (bf)';

days(2).Z     = Z_af;
days(2).t     = t_af;
days(2).label = 'Day 2 (af)';

if ~isempty(Z_day3) && ~isempty(t_day3)
    days(3).Z     = Z_day3;
    days(3).t     = t_day3;
    days(3).label = 'Day 3 (fasted)';
end

nDays = numel(days);

% ---------- compute a common cell sorting based on odor 2, day 1 ----------
ord = [];
if numel(Z_bf) >= 2 && ~isempty(Z_bf{2})
    Zref = Z_bf{2};                         % odor 2, day 1
    Mref = mean(Zref, 3, 'omitnan');        % time x cells
    win  = BS*FR+1 : (BS+2)*FR;             % e.g. 3–5 s after onset
    resp_mean = mean(Mref(win, :), 1, 'omitnan');
    [~, ord] = sort(resp_mean, 'descend');  % strongest responders to odor2 on top
else
    warning('plot_heatmaps_multiday: Odor 2 on day 1 is missing/empty. Using unsorted order.');
end

% ---------- maximum number of odors across days ----------
K_bf   = numel(Z_bf);
K_af   = numel(Z_af);
K_day3 = numel(Z_day3);
K = max([K_bf, K_af, K_day3]);

% ---------- plot for each odor k ----------
for k = 1:K
    % skip if this odor doesn't exist in ANY day
    has_any = false;
    for d = 1:nDays
        if k <= numel(days(d).Z) && ~isempty(days(d).Z{k})
            has_any = true;
            break;
        end
    end
    if ~has_any, continue; end

    fig = figure('Color','w', 'Name', sprintf('Odor %d (%s)', k, tag), ...
                 'Position',[100 100 600 200*nDays]);

    for d = 1:nDays
        if k > numel(days(d).Z) || isempty(days(d).Z{k})
            continue;
        end

        Z = days(d).Z{k};        % time x cells x trials
        t = days(d).t(:);        % time vector (column)

        M = mean(Z, 3, 'omitnan');   % time x cells

        % apply common sorting from odor2 day1 if available
        if ~isempty(ord) && size(M,2) >= max(ord)
            M = M(:, ord);
        end

        subplot(nDays, 1, d);
        h = heatmap(M', 'Colormap', parula, 'GridVisible', 'off');
        caxis([-1 2]);

        % Add vertical odor on/off lines using the time vector of that day
        S  = struct(h);
        ax = S.Axes;

        % find indices closest to 0 s and 2 s
        [~, idx_on]  = min(abs(t - 0));
        [~, idx_off] = min(abs(t - 2));
        xline(ax, idx_on + 0.5,  'k', 'LineWidth', 2);
        xline(ax, idx_off + 0.5, 'k', 'LineWidth', 2);

        title(sprintf('%s — Odor %d', days(d).label, k));
    end

    sgtitle(sprintf('Heatmaps (sorted by odor 2, Day 1) — %s — Odor %d', tag, k));

    fname = sprintf('heatmap_multiday_%s_odor%02d.png', tag, k);
    saveas(fig, fullfile(figdir, fname));
    close(fig);
end
end
