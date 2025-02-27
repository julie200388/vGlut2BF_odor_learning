dprime_mouse_updated_all=zeros(10,n);
for mouse=1:n
    if mouse==6
        continue
    end
   
kall=[];
for i=1:size(mousenum{1,mouse}.data,2)
    k=convertCharsToStrings(mousenum{1,mouse}.data(i).subject);
    kall=cat(1,kall,k);
end

% Convert kall to categorical if needed for unique processing
[unique_vals, ~, ic] = unique(kall); 

% Extract dprime values
dprime_mouse = zeros(length(kall),1);
for i=1:size(mousenum{1,mouse}.data,2)
    dprime_mouse(i) = mousenum{1,mouse}.data(i).dprime;
end
dprime_mouse(isnan(dprime_mouse)) = 0;  
% Identify unique values that appear more than once
[unique_counts, ~] = histcounts(ic, 1:max(ic)+1);
repeated_values = unique_vals(unique_counts > 1);

% Compute mean dprime for each repeated value and update dprime_mouse
for i = 1:length(repeated_values)
    idx = find(kall == repeated_values(i)); % Find all occurrences
    mean_dprime = nanmean(dprime_mouse(idx)); % Compute mean
    dprime_mouse(idx) = mean_dprime; % Assign mean to all occurrences
end

% Remove duplicate dprime values while maintaining order
dprime_mouse_updated = unique(dprime_mouse,'stable');
dprime_mouse_updated_all(:,mouse)=dprime_mouse_updated;
clear dprime_mouse_updated
end