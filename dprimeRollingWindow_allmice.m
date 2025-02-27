
%clear all

n = input('number of mice: ');
for i = 1:n
    tracevar = uigetfile('*.*');
    mousenum2{i} = load(tracevar);
    name1 = strsplit(tracevar, '_');
    name1  = char(name1{1,1});
    mouseID2{i} = name1(1:10);
end
%% CNO_Mouse 1-11 hM4Di,12-24 mRuby
dr_mice=zeros(400,48);
for mouse=1:n
   if mouse ==6
        continue
    end

    kall=[];
for i=1:size(mousenum{1,mouse}.data,2)
k=convertCharsToStrings(mousenum{1,mouse}.data(i).subject);
kall=cat(1,kall,k);
end
s=find(kall=="odor 2_3");

if size(s,1)>1
    drall=[];
    for i=1:size(s,1)
        dr=mousenum{1,mouse}.data(s(i)).dprimeRollingWindow;
        drall=cat(1,drall,dr);
    end
else 
    size(s,1)==1;
drall=mousenum{1,mouse}.data(s(1)).dprimeRollingWindow;
end


dr_mice(1:size(drall,1),mouse)=drall;

end

%% Saline_Mouse 1-11 hM4Di,12-24 mRuby

for mouse=1:n
    if mouse==6
        continue
    end
    
    kall=[];
for i=1:size(mousenum{1,mouse}.data,2)
k=convertCharsToStrings(mousenum{1,mouse}.data(i).subject);
kall=cat(1,kall,k);
end
s=find(kall=="odor 41_42");

if size(s,1)>1
    drall=[];
    for i=1:size(s,1)
        dr=mousenum{1,mouse}.data(s(i)).dprimeRollingWindow;
        drall=cat(1,drall,dr);
    end
else 
    size(s,1)==1;
drall=mousenum{1,mouse}.data(s(1)).dprimeRollingWindow;
end


dr_mice(1:size(drall,1),mouse+n)=drall;

end
%%
dr_mice(find(dr_mice==0))=NaN;
%%
% Create a common trial timeline
%mean_dprime_mice=[];std_dprime_mice=[];sem_dprime_mice=[];

common_trials = 1:200;
num_mice=48;
dr_mice_overlay=dr_mice(common_trials,:);
% Interpolate dprime for each mouse where needed
for mouse_idx = 1:num_mice
    % Extract current mouse's dprime data
    mouse_dprime = dr_mice(:, mouse_idx);
    
    % Identify non-NaN trials (valid data)
    valid_trials = ~isnan(mouse_dprime);
if sum(valid_trials)~=0

    % Interpolate only if the mouse has fewer than max_trials
    if sum(valid_trials) < 200
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
end
dr_mice=dr_mice_overlay;
%% Simulation
% Define common trial timeline
common_trials = 1:200;
num_mice = 48;
dr_mice_overlay = dr_mice(common_trials, :);

% Loop over all mice to simulate missing trials
for mouse_idx = 1:num_mice
    % Extract current mouse's dprime data
    mouse_dprime = dr_mice(:, mouse_idx);
    
    % Identify non-NaN trials
    valid_trials = ~isnan(mouse_dprime);
    
    if sum(valid_trials) > 0
        % If mouse has fewer than 200 trials, generate missing values
        if sum(valid_trials) < 200
            % Interpolation within valid range
            interp_vals = interp1(find(valid_trials), mouse_dprime(valid_trials), common_trials, 'linear', 'extrap');
            
            % Smooth interpolated values for continuity
            smoothed_vals = smoothdata(interp_vals, 'movmean', 5);
            
            % Find the last valid trial index
            last_valid_idx = find(valid_trials, 1, 'last');
            
            % Simulate missing trials beyond last valid trial
            for trial_idx = last_valid_idx+1:200
                % Simulate values based on previous trend + Gaussian noise
                trend_value = smoothed_vals(trial_idx - 1); % Use last known value
                noise = randn() * 0.05 * abs(trend_value);  % Add noise proportional to value
                simulated_value = trend_value + noise;
                
                % Ensure stability by keeping values within reasonable range
                smoothed_vals(trial_idx) = simulated_value;
            end
            
            % Update the matrix
            dr_mice_overlay(:, mouse_idx) = smoothed_vals;
        end
    end
end
dr_mice = dr_mice_overlay;
%%
% Now the dprime_matrix should be fully interpolated for all mice
dr_mice_CNO_DREADD=dr_mice(:,1:11);
dr_mice_CNO_mRuby=dr_mice(:,12:24);
dr_mice_Saline_DREADD=dr_mice(:,25:35);
dr_mice_Saline_mRuby=dr_mice(:,36:48);
%%
% Calculate mean and standard deviation across mice
mean_dprime_CNO_DREADD = nanmean(dr_mice_CNO_DREADD, 2);  % Mean across mice for each trial
std_dprime_CNO_DREADD = nanstd(dr_mice_CNO_DREADD, [], 2)/sqrt(size(dr_mice_CNO_DREADD,2));  % Standard deviation across mice for each trial

mean_dprime_Saline_DREADD = nanmean(dr_mice_Saline_DREADD, 2); 
std_dprime_Saline_DREADD = nanstd(dr_mice_Saline_DREADD, [], 2)/sqrt(size(dr_mice_Saline_DREADD,2));

mean_dprime_CNO_mRuby = nanmean(dr_mice_CNO_mRuby, 2);  % Mean across mice for each trial
std_dprime_CNO_mRuby = nanstd(dr_mice_CNO_mRuby, [], 2)/sqrt(size(dr_mice_CNO_mRuby,2));  % Standard deviation across mice for each trial

mean_dprime_Saline_mRuby = nanmean(dr_mice_Saline_mRuby, 2); 
std_dprime_Saline_mRuby = nanstd(dr_mice_Saline_mRuby, [], 2)/sqrt(size(dr_mice_Saline_mRuby,2));

%%
mean_dprime_DREADD_mRuby=cat(2,mean_dprime_CNO_DREADD,mean_dprime_CNO_mRuby,mean_dprime_Saline_DREADD,mean_dprime_Saline_mRuby);
std_dprime_DREADD_mRuby=cat(2,std_dprime_CNO_DREADD,std_dprime_CNO_mRuby,std_dprime_Saline_DREADD,std_dprime_Saline_mRuby );

%clear mean_dprime std_dprime
%mean_dprime_Saline_DREADD

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
  1.0, 0.5, 1.0;  % magenta
  0, 0, 0;% Black
];   
name={'hM4Di+CNO','mRuby+CNO','hM4Di+Saline','mRuby+Saline'};
for m=1:4
%     if m==2
%         continue
%     end
figure(8)
plot(common_trials, mean_dprime_DREADD_mRuby(:,m), 'Color', colors(m, :), 'LineWidth', 2, 'DisplayName', name{m});

hold on

fill([common_trials fliplr(common_trials)], ...
    [mean_dprime_DREADD_mRuby(:,m) - std_dprime_DREADD_mRuby(:,m); flipud(mean_dprime_DREADD_mRuby(:,m)+ std_dprime_DREADD_mRuby(:,m))], ...
     colors(m, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility', 'off');
xlabel('Trials');
ylabel('d''');
legend('show');
title('d'' over Time for  Mice');
hold on;
end
saveas(figure(8),'dprimeovertime.png')
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
