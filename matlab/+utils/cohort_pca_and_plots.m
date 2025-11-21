function cohort_pca_and_plots(MeanPerCell_sets, figdir, tag, varargin)
% ... (docstring unchanged)

% --------- Parse inputs ----------
p = inputParser;
addParameter(p, 'Colors', [], @(x) isempty(x) || (isnumeric(x)&&size(x,2)==3));
addParameter(p, 'SetLabels', {});
addParameter(p, 'OdorIdx', [], @(x) isempty(x) || isnumeric(x));   % <-- NEW
parse(p, varargin{:});
clr_user = p.Results.Colors;
labels   = p.Results.SetLabels;
OdorIdx  = p.Results.OdorIdx;

S = numel(MeanPerCell_sets);
K_all = min(cellfun(@numel, MeanPerCell_sets));

% Default to all odors if none specified:
if isempty(OdorIdx), OdorIdx = 1:K_all; end
K = numel(OdorIdx);  % number of odors to plot

% Default colors
default_clr = [
    0,    1,    0;
    0.92, 0.61, 0.61;
    1.0,  0.5,  0.0;
    0.06, 1,    1;
    0.75, 0.60, 1.0;
    0.0,  0.6,  0.9;
    0.4,  0.0,  0.6
];
if isempty(clr_user)
    clr = default_clr(1:K, :);
else
    clr = clr_user(1:K, :);
end
if isempty(labels)
    labels = arrayfun(@(s) sprintf('set%d', s), 1:S, 'uni', 0);
end

% --------- Build PCA matrix using ONLY selected odors ----------
M = [];
segments = struct('set',{},'odor',{},'r1',{},'r2',{});
for s = 1:S
    for ii = 1:K
        k = OdorIdx(ii);                % map into selected odor index
        X = MeanPerCell_sets{s}{k};     % time x cells
        if isempty(X), continue; end
        T = size(X,1);
        r1 = size(M,1) + 1;
        M  = [M; X]; %#ok<AGROW>
        r2 = size(M,1);
        segments(end+1) = struct('set',s,'odor',ii,'r1',r1,'r2',r2); %#ok<AGROW>
        % NOTE: segments.odor is 1..K within the selected subset, so colors match.
    end
end
if isempty(M)
    warning('No data found for PCA plotting (OdorIdx=%s).', mat2str(OdorIdx));
    return;
end

% --------- PCA (cells as variables) ----------
CV = cov(M);
[ev, el] = eig(CV);
[~, order] = sort(diag(el), 'ascend');
ev = ev(:, order);
data_3d = M * ev(:, end-2:end);

% --------- Plot ----------
fig = figure('Color','w'); hold on;

darkenFactor = 0.6;  % 0.5–0.8; smaller = darker

for s = 1:S
    for ii = 1:K
        idx = find([segments.set]==s & [segments.odor]==ii, 1);
        if isempty(idx), continue; end
        r = segments(idx).r1:segments(idx).r2;

        segX = data_3d(r,3); segY = data_3d(r,2); segZ = data_3d(r,1);

        base_color = clr(ii,:);
        if s == 1
            line_color = base_color;                % BEFORE: original color
        else
            line_color = base_color * darkenFactor; % AFTER: darker same color
        end

        % trajectory
        plot3(segX, segY, segZ, '-', ...
              'Color', line_color, 'LineWidth', 2.5, ...
              'DisplayName', sprintf('%s: odor %d', labels{s}, OdorIdx(ii)));

        % markers
        plot3(segX(1), segY(1), segZ(1), 'o', 'MarkerSize', 8, ...
              'MarkerFaceColor', line_color, 'MarkerEdgeColor','none', 'HandleVisibility','off');
        plot3(segX(round(end/2)), segY(round(end/2)), segZ(round(end/2)), 'd', 'MarkerSize', 8, ...
              'MarkerFaceColor', line_color, 'MarkerEdgeColor','none', 'HandleVisibility','off');
        plot3(segX(end), segY(end), segZ(end), '^', 'MarkerSize', 8, ...
              'MarkerFaceColor', line_color, 'MarkerEdgeColor','none', 'HandleVisibility','off');
    end
end

xl=xlim; yl=ylim; zl=zlim;
xticks([floor(xl(1)) ceil(xl(2))]); yticks([floor(yl(1)) ceil(yl(2))]); zticks([floor(zl(1)) ceil(zl(2))]);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
view([210 30]); grid on; box off;
ax=gca; ax.LineWidth=2; ax.FontSize=14;
legend('Location','northeast','Box','off');
title(sprintf('Cohort PCA trajectories (%s) — odors %s', tag, mat2str(OdorIdx)));

saveas(fig, fullfile(figdir, sprintf('cohort_PCA_%s_odors_%s.png', tag, strrep(mat2str(OdorIdx),' ','_'))));
close(fig);
