% +utils/subtract_baseline_condition.m
function Sub = subtract_baseline_condition(Z_byOdor, resp, mo_index)
% subtract MO (odor mo_index) traces for cells flagged RC/IC separately
Sub = Z_byOdor;
if mo_index<1 || mo_index>numel(Z_byOdor) || isempty(Z_byOdor{mo_index})
    return
end
MO = Z_byOdor{mo_index}; % time x cells x trials
MOmean = mean(MO,3,'omitnan'); % time x cells
for k = 1:numel(Z_byOdor)
    if isempty(Z_byOdor{k}), continue; end
    X = Z_byOdor{k};
    % RC subtraction
    if ~isempty(resp(1).RC)
        X(:,resp(1).RC,:) = X(:,resp(1).RC,:) - MOmean(:,resp(1).RC);
    end
    % IC subtraction
    if ~isempty(resp(1).IC)
        X(:,resp(1).IC,:) = X(:,resp(1).IC,:) - MOmean(:,resp(1).IC);
    end
    Sub{k} = X;
end
end
