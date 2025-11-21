function run_cohort_singleday(cfg)
% RUN_COHORT_SINGLEDAY
% Aggregate multiple animals (single-day sessions), perform:
%  - RC/IC/NC mean±CI traces (per odor, cohort)
%  - Heatmaps of all cells (per odor)
%  - Percent neurons responding to 0..K odors (|z|>3σ)
%  - PCA trajectories across cohort (optionally skipping odor 1)
% No decoding feature files are written.

% -------- Required cfg fields --------
% cfg.FR, cfg.BS, cfg.RS
% cfg.animals : array of structs, each with either:
%    trace_csv, odors_csv, onsets_csv
% or preloaded:
%    Time, AllAC, odors_vec, onset_time_sec
%
% Optional:
% cfg.drop_first_label (true/false)  % drop first odor label (e.g., MO)
% cfg.mo_index (default 1)
% cfg.subtract_mo (default false)
% cfg.out_dir (default './outputs_cohort_singleday')
% cfg.pca_odors = [] or vector (e.g., 2:K) to select odors for PCA & plots

% --- defaults ---
assert(isfield(cfg,'FR')&&isfield(cfg,'BS')&&isfield(cfg,'RS'),'Need FR/BS/RS');

if ~isfield(cfg,'mo_index'),        cfg.mo_index=1;            end
if ~isfield(cfg,'subtract_mo'),     cfg.subtract_mo=false;     end
if ~isfield(cfg,'pca_odors'),       cfg.pca_odors=[];          end

dateTag = datestr(now,'yyyymmdd_HHMMSS');
outroot = fullfile(getfield_or(cfg,'out_dir','./outputs_cohort_singleday'), dateTag);
figdir  = fullfile(outroot,'figures');
datadir = fullfile(outroot,'data');
ensure_dir(outroot); ensure_dir(figdir); ensure_dir(datadir);

FR=cfg.FR; BS=cfg.BS; RS=cfg.RS;
tvec = (-BS : 1/FR : RS);               % plotting time vector
TF   = BS*FR + RS*FR + 1;               % frames per trial window

% ===== 1) Load each animal and build Z_byOdor per animal =====
A = cfg.animals;
NA = numel(A);

All_Z_byOdor = {};      % cell(1,K) later, after we know K
All_Z_byOdor_sub = {};
Kmax = 0;

for a = 1:NA
    % --- load data for this animal ---
    if isfield(A(a),'Time') && isfield(A(a),'AllAC')
        Time  = A(a).Time(:);
        allAC = A(a).AllAC;
    else
        TBL = readtable(A(a).trace_csv,'TextType','string');
        Time = TBL.("Time_s__CellStatus");
        allAC = TBL.("accepted");
        vi=1;  % accepted_1, accepted_2...
        while true
            cname = sprintf('accepted_%d',vi);
            if any(strcmp(TBL.Properties.VariableNames,cname))
                allAC = cat(2, allAC, TBL.(cname));
                vi=vi+1;
            else
                break
            end
        end
    end
    allAC(isnan(allAC))=0;

    if isfield(A(a),'odors_vec')
        odors = A(a).odors_vec(:)'; 
    else
        odors = readmatrix(A(a).odors_csv)'; 
    end
    if isfield(A(a),'onset_time_sec')
        onsets = A(a).onset_time_sec(:)';
    else
        onsets = readmatrix(A(a).onsets_csv)'; 
    end

    if cfg.drop_first_label, odors(1)=[]; end
    onsets =  utils.load_dio_and_onsets(cfg.animals(a).onsets_csv, cfg.dio_threshold); 
    Ka = max(odors);
    Kmax = max(Kmax, Ka);

    % --- build Z_byOdor for this animal ---
    Z_byOdor_a = cell(1,Ka);
    for tr = 1:numel(onsets)
        k = odors(tr);
        [~, idx] = min(abs(Time - onsets(tr)));
        i1 = idx - BS*FR; i2 = idx + RS*FR;
        if i1<1 || i2>numel(Time), continue; end
        temp = allAC(i1:i2,:);
        base = mean(temp(FR+1:BS*FR,:),1);
        sd   = std( temp(FR+1:BS*FR,:),0,1);
        ztemp= (temp - base) ./ max(sd, eps);
        Z_byOdor_a{k} = cat(3, Z_byOdor_a{k}, ztemp);
    end

    % --- accumulate into cohort containers (concat CELLS horizontally, TRIALS along 3rd dim) ---
    if isempty(All_Z_byOdor)
        All_Z_byOdor = cell(1,Kmax);
    end
    for k = 1:Ka
        Zk = Z_byOdor_a{k};
        if isempty(Zk), continue; end
        if isempty(All_Z_byOdor{k})
            All_Z_byOdor{k} = Zk;                      % T x C x N
        else
            % pad time if needed (align to min T)
            Tmin = min(size(All_Z_byOdor{k},1), size(Zk,1));
            All_Z_byOdor{k} = All_Z_byOdor{k}(1:Tmin,:,:);
            Zk = Zk(1:Tmin,:,:);
            % concat CELLS horizontally (different animals have different C)
            All_Z_byOdor{k} = cat(2, All_Z_byOdor{k}, Zk); % T x (Csum) x Ntrials_a
        end
    end
end

K = numel(All_Z_byOdor);
% unify time length across odors
Tmin_overall = inf;
for k=1:K
    if ~isempty(All_Z_byOdor{k})
        Tmin_overall = min(Tmin_overall, size(All_Z_byOdor{k},1));
    end
end
for k=1:K
    if ~isempty(All_Z_byOdor{k})
        All_Z_byOdor{k} = All_Z_byOdor{k}(1:Tmin_overall,:,:);
    end
end
tvec = linspace(-BS, RS, Tmin_overall).';

% ===== 2) Optional MO subtraction (per-cell) and cohort means =====
if cfg.subtract_mo && cfg.mo_index>=1 && cfg.mo_index<=K && ~isempty(All_Z_byOdor{cfg.mo_index})
    mo_mean = mean(All_Z_byOdor{cfg.mo_index}, 3, 'omitnan'); % T x C
    All_Z_byOdor_sub = All_Z_byOdor;
    for k=1:K
        if k==cfg.mo_index || isempty(All_Z_byOdor{k}), continue; end
        Z = All_Z_byOdor{k};
        for tr = 1:size(Z,3)
            All_Z_byOdor_sub{k}(:,:,tr) = Z(:,:,tr) - mo_mean;
        end
    end
else
    All_Z_byOdor_sub = All_Z_byOdor;
end

MeanPerCell = cell(1,K);
for k=1:K
    if ~isempty(All_Z_byOdor{k})
        MeanPerCell{k} = mean(All_Z_byOdor{k},3,'omitnan');
    end
end
MeanPerCell_sub = cell(1,K);
for k=1:K
    if ~isempty(All_Z_byOdor_sub{k})
        MeanPerCell_sub{k} = mean(All_Z_byOdor_sub{k},3,'omitnan');
    end
end

% ===== 3) Classify RC/IC/NC on MO-subtracted cohort means =====
resp = classify_3sigma(MeanPerCell_sub, FR, BS);

% ===== 4) RC/IC/NC mean±CI plots (cohort) =====
utils.plot_traces_with_ci(All_Z_byOdor_sub, resp, tvec, figdir, 'cohort');

% ===== 5) Heatmaps (all cells, per odor) =====
% --- Determine cell sorting order based on odor 2 responses ---
if K >= 2 && ~isempty(MeanPerCell_sub{2})
    M2 = MeanPerCell_sub{2};
    ref_resp = mean(M2(BS*FR+1:(BS+2)*FR,:), 1, 'omitnan');
    [~, sort_order] = sort(ref_resp, 'descend');
else
    warning('Odor 2 not found or empty; sorting by odor 1 instead.');
    M1 = MeanPerCell_sub{1};
    ref_resp = mean(M1(BS*FR+1:(BS+2)*FR,:), 1, 'omitnan');
    [~, sort_order] = sort(ref_resp, 'descend');
end

% --- Generate heatmaps for all odors using that order ---
for k = 1:K
    M = MeanPerCell_sub{k};
    if isempty(M), continue; end
    if size(M,2) < max(sort_order)
        warning('Odor %d matrix has fewer cells than sorting reference.', k);
        continue;
    end

    H = M(:, sort_order);  % apply same sorting
    fig = figure('Color','w');
    h = heatmap(H', 'Colormap', parula, 'GridVisible', 'off');
    caxis([-1 2]);
    S = struct(h); axh = S.Axes;
    xline(axh, [BS*FR+0.5], 'k', 'LineWidth', 2);
    xline(axh, [(BS+2)*FR+0.5], 'k', 'LineWidth', 2);
    title(sprintf('Cohort heatmap | Odor %d (sorted by odor 2)', k));
    saveas(fig, fullfile(figdir, sprintf('cohort_heatmap_sortedBy2_odor%02d.png', k)));
    close(fig);
end


% ===== 6) Responsiveness 0..K odors (excluding MO) =====
Csum = 0;
for k = 2:K  % start from odor 2 to skip MO
    if ~isempty(MeanPerCell_sub{k})
        Csum = size(MeanPerCell_sub{k}, 2);
        break;
    end
end

resp_count = zeros(1, Csum);
for c = 1:Csum
    cnt = 0;
    for k = 2:K  % exclude odor 1 (MO)
        M = MeanPerCell_sub{k};
        if isempty(M), continue; end
        base_sd = std(M(1*FR+1:3*FR, c), 0, 1, 'omitnan');
        rmean   = mean(M(BS*FR+1:(BS+2)*FR, c), 1, 'omitnan');
        if rmean > 3*base_sd || rmean < -3*base_sd
            cnt = cnt + 1;
        end
    end
    resp_count(c) = cnt;
end

% Adjust histogram bins to reflect 0..(K-1) odors (since MO excluded)
validK = K - 1;  
hcounts = histcounts(resp_count, -0.5:1:(validK + 0.5));
props   = 100 * hcounts / sum(hcounts);

Ttab = table((0:validK)', hcounts(:), props(:), ...
    'VariableNames', {'n_odors','count','percent'});
writetable(Ttab, fullfile(datadir, 'cohort_responsive_counts_noMO.csv'));

fig = figure('Color','w');
bar(0:validK, props, 'FaceColor', [0.3 0.6 0.9]);
xlabel('# odors responded (|z|>3σ, excluding MO)'); 
ylabel('Percent of neurons (%)');
title('Cohort distribution of odor responsiveness (excluding MO)');
set(gca, 'LineWidth', 1.5, 'FontSize', 12); 
box off;
saveas(fig, fullfile(figdir, 'cohort_responsive_counts_noMO_bar.png'));
close(fig);


% ===== 7) PCA (cohort) using selected odors (e.g., 2:K) =====
OdorIdx = cfg.pca_odors;
if isempty(OdorIdx), OdorIdx = 1:K; end
utils.cohort_pca_and_plots({MeanPerCell_sub}, figdir, 'cohort', ...
    'SetLabels', {'single'}, 'OdorIdx', OdorIdx);

fprintf('Done. Outputs in: %s\n', outroot);
end

% ================= helpers =================
function v = getfield_or(s, fld, def)
if isfield(s,fld), v = s.(fld); else, v = def; end
end

function ensure_dir(d)
if ~exist(d,'dir'), mkdir(d); end
end

function S = classify_3sigma(MeanPerCell, FR, BS)
K = numel(MeanPerCell);
S = repmat(struct('RC',[],'IC',[],'NC',[]),1,K);
for k=1:K
    M = MeanPerCell{k}; if isempty(M), continue; end
    base = M(1*FR+1:3*FR,:);
    resp = M(BS*FR+1:(BS+2)*FR,:);
    zsd  = std(base,0,1,'omitnan');
    r    = mean(resp,1,'omitnan');
    RC = find(r >  3*zsd);
    IC = find(r < -3*zsd);
    allc = 1:size(M,2);
    NC = setdiff(allc, union(RC,IC));
    S(k).RC=RC; S(k).IC=IC; S(k).NC=NC;
end
end
