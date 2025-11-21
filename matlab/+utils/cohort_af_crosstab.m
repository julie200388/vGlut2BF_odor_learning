function OUT = cohort_af_crosstab(A, FR, BS, odors_pair)
% COHORT_AF_CROSSTAB
% Make a 3x3 cross-tab of AFTER-condition cell classes:
% Rows = status under odor A (RC/IC/NC)
% Cols = status under odor B (RC/IC/NC)
% Cells = counts (and row-wise %)
%
% A: struct array from load_animals_for_cohort (has MeanPerCell_af{odor})
% FR, BS: frame rate & baseline-sec (for the 3σ rule)
% odors_pair: [odorA odorB], e.g., [2 3] for Pentanol vs Hexanol
%
% OUT fields:
%   per_animal(k).animal_id
%   per_animal(k).counts (3x3), per_animal(k).rowperc (3x3)
%   cohort.counts (3x3), cohort.rowperc (3x3)
%   labels struct with row/col names

if nargin < 4 || isempty(odors_pair), odors_pair = [2 3]; end
oa = odors_pair(1); ob = odors_pair(2);

labels.rows = {'A_RC','A_IC','A_NC'};
labels.cols = {'B_RC','B_IC','B_NC'};

% --- classify AFTER-condition per animal using 3σ rule on mean-per-cell traces
resp_af_all = cell(numel(A),1);
for a = 1:numel(A)
    Mp = A(a).MeanPerCell_af;   % cell array per odor: time x cells
    resp_af_all{a} = local_classify_3sigma(Mp, FR, BS);
end

% --- per-animal crosstabs
PA = struct('animal_id',[],'counts',[],'rowperc',[]);
for a = 1:numel(A)
    PA(a).animal_id = string(A(a).animal_id);
    if numel(resp_af_all{a}) < max(oa,ob) || ...
       isempty(A(a).MeanPerCell_af{oa}) || isempty(A(a).MeanPerCell_af{ob})
        PA(a).counts  = zeros(3,3);
        PA(a).rowperc = zeros(3,3);
        continue
    end

    RA = resp_af_all{a}(oa);    % struct with RC, IC, NC (indices)
    RB = resp_af_all{a}(ob);

    % map A-status -> rows; B-status -> cols
    setsA = {RA.RC, RA.IC, RA.NC};
    setsB = {RB.RC, RB.IC, RB.NC};

    C = zeros(3,3);
    for r = 1:3
        for c = 1:3
            C(r,c) = numel(intersect(setsA{r}, setsB{c}));
        end
    end
    PA(a).counts = C;

    % row-wise percentages (within each A-status group)
    rowtot = sum(C,2);
    RP = zeros(3,3);
    for r = 1:3
        if rowtot(r) > 0
            RP(r,:) = 100 * C(r,:) / rowtot(r);
        end
    end
    PA(a).rowperc = RP;
end

% --- cohort aggregate (sum counts across animals)
CohC = zeros(3,3);
for a = 1:numel(PA)
    CohC = CohC + PA(a).counts;
end
rowtot = sum(CohC,2);
CohRP = zeros(3,3);
for r = 1:3
    if rowtot(r) > 0
        CohRP(r,:) = 100 * CohC(r,:) / rowtot(r);
    end
end

OUT.per_animal = PA;
OUT.cohort.counts  = CohC;
OUT.cohort.rowperc = CohRP;
OUT.labels = labels;
OUT.odors  = struct('A',oa,'B',ob);
end

% ---- local helper: same 3σ rule you use elsewhere ----
function resp = local_classify_3sigma(MeanPerCell, FR, BS)
K = numel(MeanPerCell);
resp = repmat(struct('RC',[],'IC',[],'NC',[]), 1, K);
for k = 1:K
    M = MeanPerCell{k};     % time x cells
    if isempty(M), continue; end
    base = M((1*FR+1):(3*FR), :);
    respwin = M((BS*FR+1):((BS+2)*FR), :);
    zbase = std(base, 0, 1, 'omitnan');
    rmean = mean(respwin, 1, 'omitnan');

    RC = find(rmean >  3*zbase);
    IC = find(rmean < -3*zbase);
    allcells = 1:size(M,2);
    NC = setdiff(allcells, union(RC, IC));

    resp(k).RC = RC; resp(k).IC = IC; resp(k).NC = NC;
end
end
