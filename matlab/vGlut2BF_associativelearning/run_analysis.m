function run_analysis(config)
% RUN_ANALYSIS  Publication-ready pipeline for calcium imaging associative learning
% Usage: cfg = configs.example_config(); run_analysis(cfg);

assert(isstruct(config), 'Config must be a struct. See configs/example_config.m');

rng(config.random_seed);

% ---------- OUTPUT SETUP ----------
dateTag = datestr(now,'yyyymmdd_HHMMSS');
outroot  = fullfile(config.out_dir, dateTag);
figdir   = fullfile(outroot, 'figures');
datadir  = fullfile(outroot, 'data');
utils.ensure_dir(outroot); utils.ensure_dir(figdir); utils.ensure_dir(datadir);

% ---------- LOAD TRACE TABLE ----------
tbl = utils.load_traces(config.traces_csv);
VarNames = tbl.Properties.VariableNames;

% Accepted cells matrix (cells in columns)
acceptedVar = "accepted";
assert(any(strcmp(VarNames, acceptedVar)), 'Missing column: accepted');
allAC = tbl.(acceptedVar);

% Concatenate the rest of accepted_* columns if present
accCols = startsWith(VarNames, "accepted_");
if any(accCols)
    accNames = VarNames(accCols);
    accMat = [];
    for i = 1:numel(accNames)
        accMat = cat(2, accMat, tbl.(accNames{i}));
    end
    allAC = cat(2, allAC, accMat);
end
Time = tbl.("Time_s__CellStatus");

% Fill NaNs with zero consistently
allAC = fillmissing(allAC,'constant',0);

% ---------- SPLIT SESSIONS (bf/af) BY TIME WRAP ----------
D2start = utils.split_sessions_by_time(Time, config.session_gap_threshold);
Time_bf = Time(1:D2start);
AC_bf   = allAC(1:D2start, :);
Time_af = Time(D2start+1:end);
AC_af   = allAC(D2start+1:end, :);

% ---------- DIO / ODOR ONSETS ----------
odor_onsets_bf = utils.load_dio_and_onsets(config.dio_bf_path, config.dio_threshold);
odor_onsets_af = utils.load_dio_and_onsets(config.dio_af_path, config.dio_threshold);

% ---------- LOAD ODOR LABELS (flexible) ----------
odors_bf = readmatrix(config.odors_bf_csv);
odors_af = readmatrix(config.odors_af_csv);
odors_af = odors_af(:);
if config.drop_first_label_bf, odors_bf(1)=[]; end
if config.drop_first_label_af, odors_af(1)=[]; end

% Map arbitrary labels to 1..K consistently per session
[odors_bf_mapped, bf_map] = utils.map_labels(odors_bf);
[odors_af_mapped, af_map] = utils.map_labels(odors_af);

% ---------- BUILD TRIAL TENSORS ----------
[Trials_bf, tvec_bf] = utils.build_trials(Time_bf, AC_bf, odor_onsets_bf, ...
    config.FR, config.BS, config.RS);
[Trials_af, tvec_af] = utils.build_trials(Time_af, AC_af, odor_onsets_af, ...
    config.FR, config.BS, config.RS);

% ---------- Z-SCORE PER TRIAL (baseline window) ----------
ZTrials_bf = utils.zscore_trials(Trials_bf, config.FR, config.BS);
ZTrials_af = utils.zscore_trials(Trials_af, config.FR, config.BS);

% ---------- GROUP BY ODOR ----------
K_bf = max(odors_bf_mapped);
K_af = max(odors_af_mapped);
Z_byOdor_bf = cell(1, K_bf);
Z_byOdor_af = cell(1, K_af);
for i = 1:numel(odors_bf_mapped)
    k = odors_bf_mapped(i);
    Z_byOdor_bf{k} = cat(3, Z_byOdor_bf{k}, ZTrials_bf(:,:,i));
end
for i = 1:numel(odors_af_mapped)
    k = odors_af_mapped(i);
    Z_byOdor_af{k} = cat(3, Z_byOdor_af{k}, ZTrials_af(:,:,i));
end

% ---------- MEAN TRACE PER CELL (bf/af) ----------
MeanPerCell_bf = utils.mean_traces_per_cell(Z_byOdor_bf);
MeanPerCell_af = utils.mean_traces_per_cell(Z_byOdor_af);

% ---------- CLASSIFY CELLS (3*sigma, MO channel = config.mo_index) ----------
resp_bf = utils.classify_cells_3sigma(MeanPerCell_bf, config.FR, config.BS);
resp_af = utils.classify_cells_3sigma(MeanPerCell_af, config.FR, config.BS);

%% ---------- OPTIONAL: SUBTRACT MO (by RC/IC sets from MO) ----------
if isfield(config,'do_mo_subtraction') && config.do_mo_subtraction ...
        && ~isempty(config.mo_index) && config.mo_index <= numel(MeanPerCell_bf)

    Sub_bf = utils.subtract_baseline_condition(Z_byOdor_bf, resp_bf, config.mo_index);
    Sub_af = utils.subtract_baseline_condition(Z_byOdor_af, resp_af, config.mo_index);

    % --- NEW: recompute means after subtraction ---
    MeanPerCell_sub_bf = utils.mean_traces_per_cell(Sub_bf);
    MeanPerCell_sub_af = utils.mean_traces_per_cell(Sub_af);

    % --- NEW: reclassify RC/NC/IC on subtracted data (3σ) ---
    if isfield(config,'reclassify_after_mo') && config.reclassify_after_mo
        resp_sub_bf = utils.classify_cells_3sigma(MeanPerCell_sub_bf, config.FR, config.BS);
        resp_sub_af = utils.classify_cells_3sigma(MeanPerCell_sub_af, config.FR, config.BS);

        % Save indices for record/repro
        save(fullfile(datadir,'resp_indices_after_mo.mat'), 'resp_sub_bf','resp_sub_af');
        utils.save_resp_indices_csv(resp_sub_bf, fullfile(datadir,'resp_after_mo_bf.csv'));
        utils.save_resp_indices_csv(resp_sub_af, fullfile(datadir,'resp_after_mo_af.csv'));

        % Optional plotting using the recomputed sets
        if isfield(config,'plot_after_mo_sets') && config.plot_after_mo_sets
            % utils.plot_traces_with_ci(Sub_bf, resp_sub_bf, tvec_bf, figdir, 'bf_afterMO_RCIC');
            % utils.plot_traces_with_ci(Sub_af, resp_sub_af, tvec_af, figdir, 'af_afterMO_RCIC');

            % You may also want the overlay (all cells) after MO:
            utils.plot_allcells_overlay(MeanPerCell_sub_bf, MeanPerCell_sub_af, tvec_bf, tvec_af, figdir);
        end
    end

else
  
    % No MO subtraction (or not configured)
    Sub_bf = Z_byOdor_bf; 
    MeanPerCell_sub_bf=MeanPerCell_bf;
    resp_sub_bf=resp_bf;
    Sub_af = Z_byOdor_af;
    MeanPerCell_sub_af=MeanPerCell_af ;
    resp_sub_af=resp_af;
end

%%
% ---------- EXPORT FEATURES FOR DECODER (example: 2 s after onset) ----------
feat_bf = utils.export_features_for_decoder(Sub_bf, config.decoder_window, config.FR);
writematrix(feat_bf, fullfile(datadir, 'features_bf.csv'));
feat_af = utils.export_features_for_decoder(Sub_af, config.decoder_window, config.FR);
writematrix(feat_af, fullfile(datadir, 'features_af.csv'));
% ---------- PLOTTING: TRACES + CI ----------
utils.plot_traces_with_ci(Sub_bf, resp_sub_bf, tvec_bf, figdir, 'bf');
utils.plot_traces_with_ci(Sub_af, resp_sub_af, tvec_af, figdir, 'af');
%%
% Overlay: average across all neurons (bf vs af) per odor

utils.plot_allcells_overlay(MeanPerCell_sub_bf, MeanPerCell_sub_af, tvec_bf, tvec_af, figdir);
% ---------- HEATMAPS ----------
utils.plot_heatmaps_multiday(Sub_bf, tvec_bf, ...
                             Sub_af, tvec_af, ...
                             config.FR, config.BS, ...
                             figdir, 'bf_af');


% ---------- PCA VIS (bf, af, both) ----------
utils.pca_and_plots( MeanPerCell_sub_bf(2:end) , figdir, 'bf');
utils.pca_and_plots( MeanPerCell_sub_af(2:end) , figdir, 'af', true);
utils.pca_and_plots([MeanPerCell_sub_bf(2:end), MeanPerCell_sub_af(2:end)], figdir, 'bfaf',false,true);

% ---------- DISTANCES & PLOTS ----------
utils.distances_and_plots(Sub_bf, Sub_af, tvec_bf, tvec_af, figdir);

% ---------- RELIABILITY (split-half, fixed RNG) ----------
rel_bf = utils.reliability_split_half(Sub_bf, config.num_splits, config.random_seed);
rel_af = utils.reliability_split_half(Sub_af, config.num_splits, config.random_seed);
save(fullfile(datadir,'reliability.mat'), 'rel_bf', 'rel_af');

% ---------- SAVE SUMMARY ----------
save(fullfile(datadir,'summary.mat'), ...
    'config','resp_sub_bf','resp_sub_af','MeanPerCell_sub_bf','MeanPerCell_sub_af', ...
    'Sub_bf','Sub_af','tvec_bf','tvec_af');
%% ---------- SAVE ANIMAL PACK FOR COHORT ----------
% Use whatever ID you like (e.g., filename stem)
[~, animal_id] = fileparts(config.traces_csv);

animal_pack_path = fullfile(datadir, sprintf('%s_animal_for_cohort.mat', animal_id));
save(animal_pack_path, ...
    'animal_id', ...
    'MeanPerCell_sub_bf', 'MeanPerCell_sub_af', ...       % time x cells, per odor (averaged across trials)
    'Sub_bf', 'Sub_af', ...             % time x cells x trials
    'tvec_bf', 'tvec_af');                        % time vectors

fprintf('Saved cohort pack: %s\n', animal_pack_path);

disp('✓ Analysis finished successfully.');
end
