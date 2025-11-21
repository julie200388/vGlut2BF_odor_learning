% +utils/map_labels.m
function [mapped, map] = map_labels(labels)
u = unique(labels(:)');
map = containers.Map(num2cell(u), num2cell(1:numel(u)));
mapped = arrayfun(@(x) map(x), labels);
end
