function A = load_animals_for_cohort(mat_files)
% Returns struct array A with fields:
%   animal_id, MeanPerCell_bf, MeanPerCell_af, Z_byOdor_bf, Z_byOdor_af, tvec_bf, tvec_af

A = struct('animal_id',{},'MeanPerCell_bf',{},'MeanPerCell_af',{}, ...
           'Z_byOdor_bf',{},'Z_byOdor_af',{},'tvec_bf',{},'tvec_af',{});
for i = 1:numel(mat_files)
    S = load(mat_files{i});
    A(i).animal_id       = S.animal_id;
    A(i).MeanPerCell_bf  = S.MeanPerCell_sub_bf;
    A(i).MeanPerCell_af  = S.MeanPerCell_sub_af;
    A(i).Z_byOdor_bf     = S.Sub_bf;
    A(i).Z_byOdor_af     = S.Sub_af;
    A(i).tvec_bf         = S.tvec_bf;
    A(i).tvec_af         = S.tvec_af;
end
end
