%% create a matrix with 4 conditions
dr_mice=cat(2,dr_mice_22_87_h(1:150,:),dr_mice_33_26_h(1:150,:),dr_mice_22_87(1:150,:),dr_mice_33_26(1:150,:));

%%
% Create a common trial timeline
mean_dprime_mice=[];std_dprime_mice=[];sem_dprime_mice=[];
for m =1:4
common_trials = 1:150;
num_mice=5;
dr_mice_overlay=dr_mice(:,5*(m-1)+1:5*m);
% Interpolate dprime for each mouse where needed
for mouse_idx = 1:num_mice
    % Extract current mouse's dprime data
    mouse_dprime = dr_mice_overlay(:, mouse_idx);
    
    % Identify non-NaN trials (valid data)
    valid_trials = ~isnan(mouse_dprime);
    
    % Interpolate only if the mouse has fewer than max_trials
    if sum(valid_trials) < 150
         % Interpolate inside the valid range, with linear interpolation
        interp_vals = interp1(find(valid_trials), mouse_dprime(valid_trials), common_trials, 'linear', 'extrap');
        
        % Apply smoothing to the entire interpolated curve
        smoothed_vals = smoothdata(interp_vals, 'movmean', 5); % 'movmean' with a window of 10
        
        % Update the matrix with smoothed interpolated values
        dr_mice_overlay(:, mouse_idx) = smoothed_vals;
        
        % Optionally cap the extrapolated values to avoid steep trends
        last_valid_idx = find(valid_trials, 1, 'last');
        dr_mice_overlay(last_valid_idx+1:end, mouse_idx) = smoothed_vals(last_valid_idx);  % Cap to the last valid value
    end
end

% Now the dprime_matrix should be fully interpolated for all mice

% Calculate mean and standard deviation across mice
mean_dprime = nanmean( dr_mice_overlay, 2);  % Mean across mice for each trial
std_dprime = nanstd( dr_mice_overlay, [], 2)/sqrt(num_mice);  % Standard deviation across mice for each trial
mean_dprime_mice=cat(2,mean_dprime_mice,mean_dprime);
std_dprime_mice=cat(2,std_dprime_mice,std_dprime);
clear mean_dprime std_dprime

end
%%
% Plot individual mice d' curves
figure;
hold on;
for mouse_idx = 1:num_mice
    plot(common_trials, dr_mice_overlay(:, mouse_idx), '--', 'DisplayName', ['Mouse ' num2str(mouse_idx)]);
end
%%
% % Plot the mean d' with standard deviation shading
% Define a color matrix: Gray for some mice, Red for others
colors = [
  1.0, 0.0, 0.0;  % Red
  0.5, 0.5, 0.5;
  1.0, 0.0, 0.0;
  0.5, 0.5, 0.5;% Gray
];   
for m=3:4
figure(4)
plot(common_trials, mean_dprime_mice(:,m), 'Color', colors(m, :), 'LineWidth', 2, 'DisplayName', 'Mean d''');
hold on
fill([common_trials fliplr(common_trials)], ...
    [mean_dprime_mice(:,m) - std_dprime_mice(:,m); flipud(mean_dprime_mice(:,m)+ std_dprime_mice(:,m))], ...
     colors(m, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Std dev');
xlabel('Trials');
ylabel('d''');
%legend('show');
title('d'' over Time for 5 Mice');
hold on;
end
%%
for m=1:2
figure(3)
plot(common_trials, mean_dprime_mice(:,m), 'Color', colors(m, :), 'LineWidth', 2, 'DisplayName', 'Mean d''');
hold on
fill([common_trials fliplr(common_trials)], ...
    [mean_dprime_mice(:,m) - std_dprime_mice(:,m); flipud(mean_dprime_mice(:,m)+ std_dprime_mice(:,m))], ...
     colors(m, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Std dev');
xlabel('Trials');
ylabel('d''');
%legend('show');
title('d'' over Time for 5 Mice');
hold on;
end
%%
saveas(figure(3),'hM4di_Saline_CNO.png')
saveas(figure(4),'mRuby_Saline_CNO.png')
%%
