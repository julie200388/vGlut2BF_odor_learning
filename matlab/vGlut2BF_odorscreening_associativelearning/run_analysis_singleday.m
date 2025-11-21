function run_analysis_singleday(config)
% RUN_ANALYSIS_SINGLEDAY  One-day analysis: RC/IC/NC means, heatmaps, decoding file,
%                         and odor-responsivity distribution (0..k odors).
%
% Required config fields:
%   config.FR                 % frames/sec, e.g., 15
%   config.BS                 % baseline seconds before odor, e.g., 3
%   config.RS                 % response seconds after onset (for plotting windows), e.g., 10
%   config.drop_first_label   % true/false: drop first odor label (e.g., MO) from analyses
%
%   Input files OR preloaded tables:
%   EITHER provide a CSV path with traces in 'accepted' columns + 'Time_s__CellStatus':
%     config.trace_csv        % path to the per-frame trace table (like your miniscope export)
%   OR provide directly:
%     config.Time             % Tx1 time vector
%     config.AllAC            % T x Ncells numeric matrix (NaN allowed; will be zeroed)
%
%   Event labels (odors):
%     config.odors            % 1 x Ntrials integer labels (1..K)
%     config.onset_times      % 1 x Ntrials vector of onset times (seconds, same reference as Time)
%
% Optional:
%   config.mo_index           % which odor index is MO (default 1)
%   config.subtract_mo        % true/false; subtract MO mean trace per-cell from all odors (default false)
%   config.feature_win        % [startFrame endFrame] for decoding (default [BS*FR+1, (BS+2)*FR])
%   config.out_dir            % base output folder (default './outputs_singleday/yyyymmdd_HHMMSS')
%
% Outputs:
%   Saves figures (PNG) and data CSV/MAT into out_dir.

%% ---------- I/O & Basics ----------
assert(isfield(config,'FR') && isfield(config,'BS') && isfield(config,'RS'), ...
    'Provide config.FR, config.BS, config.RS');

if ~isfield(config,'drop_first_label'), config.drop_first_label = false; end
if ~isfield(config,'mo_index'),        config.mo_index = 1;              end
if ~isfield(config,'subtract_mo'),     config.subtract_mo = false;       end

dateTag = datestr(now,'yyyymmdd_HHMMSS');
if ~isfield(config,'out_dir') || isempty(config.out_dir)
    outroot = fullfile('./outputs_singleday', dateTag);
else
    outroot = fullfile(config.out_dir, dateTag);
end
figdir  = fullfile(outroot, 'figures');
datadir = fullfile(outroot, 'data');
ensure_dir(outroot); ensure_dir(figdir); ensure_dir(datadir);

FR = config.FR; BS = config.BS; RS = config.RS;

% Load traces if only CSV provided
if ~isfield(config,'Time') || ~isfield(config,'AllAC')
    assert(isfield(config,'trace_csv'), 'Provide config.trace_csv OR Time+AllAC.');
    tbl = readtable(config.trace_csv, 'TextType','string');
    Time = tbl.("Time_s__CellStatus");
    allAC = tbl.("accepted");
    % add additional accepted_i columns if present
    vi = 1;
    while true
        col = sprintf('accepted_%d', vi);
        if any(strcmp(tbl.Properties.VariableNames, col))
            allAC = cat(2, allAC, tbl.(col));
            vi = vi+1;
        else
            break
        end
    end
else
    Time  = config.Time(:);
    allAC = config.AllAC;
end
allAC(isnan(allAC)) = 0;

assert(isfield(config,'odors') && isfield(config,'onset_times'), ...
    'Provide config.odors (1xNtrials) and config.onset_times (1xNtrials).');

odors = config.odors(:);               % 1 x Ntrials
if config.drop_first_label, odors(1)=[]; end
onsets =  utils.load_dio_and_onsets(config.dio_path, config.dio_threshold);    % 1 x Ntrials


K = max(odors);                          % # odors (after any drop)
TF = BS*FR + RS*FR + 1;                  % frames in window [ -BS .. +RS ]
twin = (-BS : 1/FR : RS);                % time vector for window (may be TF±1 depending on FR)

%%---------- Build trial tensors: Z_byOdor{k} = [T x C x Ntrials_k] ----------
Z_byOdor = cell(1,K);
for tr = 1:numel(onsets)
    k = odors(tr);
    [~, idx] = min(abs(Time - onsets(tr)));
    i1 = idx - BS*FR;
    i2 = idx + RS*FR;
    if i1 < 1 || i2 > numel(Time), continue; end
    temp = allAC(i1:i2, :);
    % z-score using baseline section (FR+1 : BS*FR) like your original code
    base    = mean(temp(FR+1:BS*FR, :), 1);
    base_sd = std( temp(FR+1:BS*FR, :), 0, 1);
    ztemp   = (temp - base) ./ max(base_sd, eps);
    Z_byOdor{k} = cat(3, Z_byOdor{k}, ztemp);
end

% Per-odor mean (time x cells)
MeanPerCell = cell(1,K);
for k = 1:K
    if ~isempty(Z_byOdor{k})
        MeanPerCell{k} = mean(Z_byOdor{k}, 3, 'omitnan');
    else
        MeanPerCell{k} = [];
    end
end

%%---------- Optional: MO subtraction on the single-day data ----------
if config.subtract_mo && config.mo_index >= 1 && config.mo_index <= K && ~isempty(MeanPerCell{config.mo_index})
    mo_mean = MeanPerCell{config.mo_index};             % time x cells
    Z_byOdor_sub = Z_byOdor;
    for k = 1:K
        if k == config.mo_index || isempty(Z_byOdor{k}), continue; end
        Z = Z_byOdor{k};
        for tr = 1:size(Z,3)
            Z_byOdor_sub{k}(:,:,tr) = Z(:,:,tr) - mo_mean;
        end
    end
else
    Z_byOdor_sub = Z_byOdor;
end

% Recompute means after (optional) subtraction
MeanPerCell_sub = cell(1,K);
for k = 1:K
    if ~isempty(Z_byOdor_sub{k})
        MeanPerCell_sub{k} = mean(Z_byOdor_sub{k}, 3, 'omitnan');
    else
        MeanPerCell_sub{k} = [];
    end
end

%% ---------- RC / IC / NC classification (3σ) per odor on MeanPerCell_sub ----------
resp = classify_3sigma(MeanPerCell_sub, FR, BS);

%% ---------- Plot mean ± CI for RC/IC/NC (per odor) ----------
% --- RC/IC/NC mean±CI plots (single day; uses MO-subtracted data if enabled) ---
utils.plot_traces_with_ci(Z_byOdor_sub, resp, twin(:), figdir, 'single');

%% ---------- Heatmap of all cells (per odor, ordered by response magnitude) ----------
% --- Determine sorting order based on the 2nd odor ---
if K >= 2 && ~isempty(MeanPerCell_sub{2})
    Mref = MeanPerCell_sub{2};
    ref_resp = mean(Mref(BS*FR+1:(BS+2)*FR,:), 1, 'omitnan');
    [~, order_ref] = sort(ref_resp, 'descend');
else
    warning('Odor 2 not found or empty; sorting by odor 1 instead.');
    Mref = MeanPerCell_sub{1};
    ref_resp = mean(Mref(BS*FR+1:(BS+2)*FR,:), 1, 'omitnan');
    [~, order_ref] = sort(ref_resp, 'descend');
end

% --- Generate heatmaps for all odors using that order ---
for k = 1:K
    M = MeanPerCell_sub{k};
    if isempty(M), continue; end
    
    % Apply same order determined from odor 2
    if size(M,2) >= numel(order_ref)
        H = M(:, order_ref);
    else
        warning('Odor %d matrix has fewer cells than reference sorting.', k);
        H = M;
    end

    % --- Plot heatmap ---
    fig = figure('Color','w');
    h = heatmap(H', 'Colormap', parula, 'GridVisible', 'off');
    caxis([-1 2]);

    % Add odor on/off vertical lines
    S = struct(h); axh = S.Axes;
    xline(axh, [BS*FR+0.5], 'k', 'LineWidth', 2);
    xline(axh, [(BS+2)*FR+0.5], 'k', 'LineWidth', 2);

    % Save
    title(sprintf('Heatmap (all cells, sorted by odor 2) | Odor %d', k));
    saveas(fig, fullfile(figdir, sprintf('single_heatmap_sortedBy2_odor%02d.png', k)));
    close(fig);
end

%% ---------- Decoding feature file (flattened window) ----------
if ~isfield(config,'feature_win') || isempty(config.feature_win)
    feature_win = [BS*FR+1, (BS+2)*FR];
else
    feature_win = config.feature_win;
end

X = []; y = [];   % rows = trials, columns = flattened time x cells
for k = 1:K
    Z = Z_byOdor_sub{k}; if isempty(Z), continue; end
    T = feature_win(1):feature_win(2);
    for tr = 1:size(Z,3)
        feat = reshape(Z(T,:,tr), 1, []);
        X = [X; feat]; %#ok<AGROW>
        y = [y; k];    %#ok<AGROW>
    end
end

feat_path = fullfile(datadir, 'single_decoding_features.csv');
writematrix([X y], feat_path);
fprintf('Wrote decoding features: %s  (size %d x %d)\n', feat_path, size([X y],1), size([X y],2));

%% ---------- % of neurons responding to 0..(K-1) odors (excluding MO) ----------
% Decide which odor index is MO
if exist('config','var') && isfield(config,'mo_index') && ~isempty(config.mo_index)
    mo_idx = config.mo_index;    % e.g., 1
else
    mo_idx = 1;                  % default to 1 if not set
end

% Odors to evaluate (skip MO)
odor_idx = setdiff(1:K, mo_idx);
validK   = numel(odor_idx);      % number of non-MO odors

C = size(allAC, 2);
resp_count = zeros(1, C);        % how many non-MO odors each neuron responded to (RC or IC)

for c = 1:C
    count = 0;
    for k = odor_idx
        if isempty(MeanPerCell_sub{k}), continue; end
        M = MeanPerCell_sub{k};
        base_sd = std(M(1*FR+1:3*FR, c), 0, 1, 'omitnan');
        rmean   = mean(M(BS*FR+1:(BS+2)*FR, c), 1, 'omitnan');
        if rmean > 3*base_sd || rmean < -3*base_sd
            count = count + 1;
        end
    end
    resp_count(c) = count;
end

% histogram (0..validK)
hcounts = histcounts(resp_count, -0.5:1:(validK + 0.5));
props   = 100 * (hcounts / max(1, sum(hcounts)));

T = table((0:validK)', hcounts(:), props(:), ...
    'VariableNames', {'n_odors','count','percent'});
writetable(T, fullfile(datadir, 'single_responsive_counts_noMO.csv'));

fig = figure('Color','w');
bar(0:validK, props, 'FaceColor', [0.3 0.6 0.9]);
xlabel(sprintf('# odors responded (|z|>3σ, excluding odor %d as MO)', mo_idx));
ylabel('Percent of neurons (%)');
title('Distribution of odor responsiveness per neuron (single day, no MO)');
set(gca, 'LineWidth', 1.5, 'FontSize', 12); box off
saveas(fig, fullfile(figdir, 'single_responsive_counts_noMO_bar.png'));
close(fig);

fprintf('Done. Outputs in: %s (no MO in responsiveness)\n', outroot);


end % ----------------------- end main -----------------------

%% ========== Helpers ==========

function ensure_dir(d)
if ~exist(d,'dir'), mkdir(d); end
end

function S = classify_3sigma(MeanPerCell, FR, BS)
K = numel(MeanPerCell);
S = repmat(struct('RC',[],'IC',[],'NC',[]), 1, K);
for k = 1:K
    M = MeanPerCell{k};
    if isempty(M), continue; end
    base = M(1*FR+1:3*FR, :);
    resp = M(BS*FR+1:(BS+2)*FR, :);
    zsd  = std(base, 0, 1, 'omitnan');
    r    = mean(resp, 1, 'omitnan');
    RC = find(r >  3*zsd);
    IC = find(r < -3*zsd);
    allc = 1:size(M,2);
    NC = setdiff(allc, union(RC, IC));
    S(k).RC = RC; S(k).IC = IC; S(k).NC = NC;
end
end

function [m, lo, hi] = mean_ci_95(M)
% M: time x cells subset
m  = mean(M, 2, 'omitnan');
se = std(M, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(M),2));
ci = 1.96 * se;
lo = m - ci; hi = m + ci;
end
%% Howto call it
% cfg.FR = 15;           % Hz
% cfg.BS = 3;            % sec before odor
% cfg.RS = 10;           % sec after odor
% 
% cfg.drop_first_label = true;      
% cfg.mo_index = 1;      
% cfg.subtract_mo = false;          
% 
% cfg.feature_win = [];             
% 
% cfg.trace_csv  = '/Users/peyshyuanchin/Documents/Project_AON/092525_4odor_H1_traces_28cells.csv';   % your calcium trace file
% 
% % --- Load odors and onsets from CSV files ---
% cfg.odors       = readmatrix('/Users/peyshyuanchin/Documents/Project_AON/6stim_10trials.csv');     % e.g., [1 2 3 4 5 6 1 2 3 ...]
% cfg.dio_path = '/Users/peyshyuanchin/Documents/Project_AON/092525_4odor_H1_GPIO.csv';             % e.g., [10.2 20.3 30.4 ...] in seconds
% cfg.dio_threshold       = 1000;   % threshold to binarize DIO
% 
% cfg.out_dir = './outputs_singleday';

% run_analysis_singleday(cfg);