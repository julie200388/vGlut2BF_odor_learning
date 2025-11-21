function save_resp_indices_csv(resp, outcsv)
% SAVE_RESP_INDICES_CSV  Save RC/NC/IC indices per odor to CSV
% resp: 1xK struct array with fields RC, IC, NC (row vectors of indices)

K = numel(resp);
rows = {};
for k = 1:K
    RC = sprintf('%s', strjoin(arrayfun(@num2str, resp(k).RC, 'uni',0), ' '));
    IC = sprintf('%s', strjoin(arrayfun(@num2str, resp(k).IC, 'uni',0), ' '));
    NC = sprintf('%s', strjoin(arrayfun(@num2str, resp(k).NC, 'uni',0), ' '));
    rows(end+1,:) = {k, RC, IC, NC}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Odor','RC_indices','IC_indices','NC_indices'});
writetable(T, outcsv);
end
