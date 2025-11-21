% +utils/mean_traces_per_cell.m
function MeanPerCell = mean_traces_per_cell(Z_byOdor)
MeanPerCell = cell(size(Z_byOdor));
for k = 1:numel(Z_byOdor)
    if isempty(Z_byOdor{k}), MeanPerCell{k} = []; continue; end
    MeanPerCell{k} = mean(Z_byOdor{k}, 3, 'omitnan'); % time x cells
end
end
