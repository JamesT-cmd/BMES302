% Load the table containing the calibration data
load('calibration2_data_table.mat');

% Extract the data for 'topToBottom.txt'
topToBottom_data = [];
for i = 1:height(data_table)
    if strcmp(data_table.Filename{i}, 'topToBottom.txt')
        topToBottom_data = data_table.Data{i};
        break;
    end
end

% Check if data is found
if isempty(topToBottom_data)
    error('Data for topToBottom.txt not found in data_table.');
end

% Initialize arrays to store the data
leftmost_positive_spikes = [];
rightmost_negative_spikes = [];

% Process the data using the process_data function
% Assuming process_data returns processed data
calibrated_data = process_data(topToBottom_data, 200, 2, 'topToBottom.txt');

% Calculate the mean and standard deviation of the entire calibrated data
mean_data = mean(calibrated_data);
std_data = std(calibrated_data);
threshold_positive = mean_data + std_data;
threshold_negative = mean_data - std_data; % Using negative threshold for negative peaks

% Identify the positive peaks in the entire data
[pks_positive, locs_positive] = findpeaks(calibrated_data);
significant_locs_positive = locs_positive(pks_positive > threshold_positive);

% Identify the negative peaks in the entire data
[pks_negative, locs_negative] = findpeaks(-calibrated_data);
significant_locs_negative = locs_negative(pks_negative > threshold_positive);

% Group significant locations into chunks
chunk_threshold = 50; % Define a threshold for grouping points into chunks

% Group positive spikes
if ~isempty(significant_locs_positive)
    grouped_locs_positive = {significant_locs_positive(1)};
    for i = 2:length(significant_locs_positive)
        if significant_locs_positive(i) - significant_locs_positive(i-1) <= chunk_threshold
            grouped_locs_positive{end} = [grouped_locs_positive{end}, significant_locs_positive(i)];
        else
            grouped_locs_positive{end+1} = significant_locs_positive(i);
        end
    end
end

% Group negative spikes
if ~isempty(significant_locs_negative)
    grouped_locs_negative = {significant_locs_negative(1)};
    for i = 2:length(significant_locs_negative)
        if significant_locs_negative(i) - significant_locs_negative(i-1) <= chunk_threshold
            grouped_locs_negative{end} = [grouped_locs_negative{end}, significant_locs_negative(i)];
        else
            grouped_locs_negative{end+1} = significant_locs_negative(i);
        end
    end
end

% Get the leftmost side of the positive spikes
for i = 1:length(grouped_locs_positive)
    group = grouped_locs_positive{i};
    leftmost_positive_spikes = [leftmost_positive_spikes; group(1)];
end

% Get the rightmost side of the negative spikes
for i = 1:length(grouped_locs_negative)
    group = grouped_locs_negative{i};
    rightmost_negative_spikes = [rightmost_negative_spikes; group(end)];
end

% Convert indices to time
time_leftmost_positive_spikes = leftmost_positive_spikes / 200;
time_rightmost_negative_spikes = rightmost_negative_spikes / 200;

% Plot the data and identified significant spikes for visualization
figure;
plot((1:length(calibrated_data)) / 200, calibrated_data);
hold on;

% Plot the leftmost side of positive spikes and rightmost side of negative spikes
plot(time_leftmost_positive_spikes, calibrated_data(leftmost_positive_spikes), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Left Side of Positive Peaks');
plot(time_rightmost_negative_spikes, calibrated_data(rightmost_negative_spikes), 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Right Side of Negative Peaks');


yline(26.57, 'b--', 'DisplayName', '26.57 Degrees');
yline(-26.57, 'b--', 'DisplayName', '-26.57 Degrees');
yline(40.88, 'r--', 'DisplayName', '40.88 Degrees');
yline(-40.88, 'r--', 'DisplayName', '-40.88 Degrees');
yline(45, 'g--', 'DisplayName', '45 Degrees')
yline(-45, 'g--', 'DisplayName', '-45 Degrees')

title('Processed EOG Data with Identified Significant Peaks and Edges');
xlabel('Time (s)');
ylabel('EOG Signal');
legend('EOG Data', 'Left Side of Positive Peaks', 'Right Side of Negative Peaks', '26.57 Degrees', '-26.57 Degrees', '40.88 Degrees', '-40.88 Degrees');
hold off;

% Create a table to display the identified values
T = table(time_leftmost_positive_spikes, calibrated_data(leftmost_positive_spikes), ...
          time_rightmost_negative_spikes, calibrated_data(rightmost_negative_spikes), ...
          'VariableNames', {'Time_Leftmost_Positive_Spike', 'Leftmost_Positive_Spike_Value', ...
                            'Time_Rightmost_Negative_Spike', 'Rightmost_Negative_Spike_Value'});

% Display the table
disp(T);

angle_pairs = [0, 26.5; 26.5, 41];

T = addvars(T, [41, 41, 26.5, 26.5, 0, 0, -26.5, -26.5, -41, -41]', 'NewVariableNames', 'ZY_angle')
% Define the desired order
desired_order = [5 6 3 4 7 8 1 2 9 10];

% Rearrange the rows of the table according to the desired order
T_reordered = T(desired_order, :);

% Display the reordered table
disp(T_reordered);

% Calculate the magnitudes of angles
T.ZY_angle_magnitude = abs(T.ZY_angle);

% Calculate the magnitudes of spike values
T.Leftmost_Positive_Spike_Magnitude = abs(T.Leftmost_Positive_Spike_Value);
T.Rightmost_Negative_Spike_Magnitude = abs(T.Rightmost_Negative_Spike_Value);

% Extract unique angles from the angle_pairs
unique_angles = unique(angle_pairs);

% Initialize arrays to store means, standard errors, confidence intervals, and p-values
means = zeros(length(unique_angles), 1);
standard_errors = zeros(length(unique_angles), 1);
conf_intervals = zeros(length(unique_angles), 2);
p_values_combined = zeros(size(angle_pairs, 1), 1);

alpha = 0.05; % significance level for confidence intervals
t_value = tinv(1 - alpha/2, 3); % t-value for 95% confidence interval with 3 degrees of freedom

% Calculate mean, standard error, and confidence intervals for each unique angle
for i = 1:length(unique_angles)
    current_angle = unique_angles(i);
    
    % Extract data for the current angle magnitude
    data_current_angle_pos = T.Leftmost_Positive_Spike_Magnitude(T.ZY_angle_magnitude == current_angle);
    data_current_angle_neg = T.Rightmost_Negative_Spike_Magnitude(T.ZY_angle_magnitude == current_angle);
    
    % Combine the positive and negative spike values
    data_current_angle_combined = [data_current_angle_pos; data_current_angle_neg];
    
    % Calculate mean, standard error, and confidence intervals
    means(i) = mean(data_current_angle_combined);
    standard_errors(i) = std(data_current_angle_combined) / sqrt(length(data_current_angle_combined));
    conf_intervals(i, :) = means(i) + t_value * [-standard_errors(i), standard_errors(i)];
end

% Perform statistical tests for specified angle pairs
for i = 1:size(angle_pairs, 1)
    angle1 = angle_pairs(i, 1);
    angle2 = angle_pairs(i, 2);
    
    % Extract data for the two angles
    data_angle1_pos = T.Leftmost_Positive_Spike_Magnitude(T.ZY_angle_magnitude == angle1);
    data_angle1_neg = T.Rightmost_Negative_Spike_Magnitude(T.ZY_angle_magnitude == angle1);
    data_angle2_pos = T.Leftmost_Positive_Spike_Magnitude(T.ZY_angle_magnitude == angle2);
    data_angle2_neg = T.Rightmost_Negative_Spike_Magnitude(T.ZY_angle_magnitude == angle2);
    
    % Combine the positive and negative spike values for both angles
    data_angle1_combined = [data_angle1_pos; data_angle1_neg];
    data_angle2_combined = [data_angle2_pos; data_angle2_neg];
    
    % Perform paired t-test for combined magnitudes
    if length(data_angle1_combined) >= 4 && length(data_angle2_combined) >= 4
        [~, p_combined] = ttest(data_angle1_combined(1:4), data_angle2_combined(1:4));
        p_values_combined(i) = p_combined;
    else
        p_values_combined(i) = NaN;
    end
end

% Display results
result_table = table(angle_pairs(:, 1), angle_pairs(:, 2), p_values_combined, ...
                     'VariableNames', {'Angle1', 'Angle2', 'P_Value_Combined_Magnitude'});
disp(result_table);

% Plot means with error bars and confidence intervals
figure;
hold on;

% Define colors
colors = lines(length(unique_angles));

% Plot confidence intervals
for i = 1:length(unique_angles)
    fill([i-0.2, i+0.2, i+0.2, i-0.2], [conf_intervals(i, 1), conf_intervals(i, 1), conf_intervals(i, 2), conf_intervals(i, 2)], ...
         colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Plot means with error bars
for i = 1:length(unique_angles)
    errorbar(i, means(i), standard_errors(i), 'o', 'Color', colors(i,:), 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', sprintf('Angle Magnitude %.1f', unique_angles(i)));
    text(i, means(i) + standard_errors(i) + 1, sprintf('%.2f', means(i)), 'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', colors(i,:));
end

% Add significance stars and lines for specified angle pairs
for i = 1:size(angle_pairs, 1)
    if ~isnan(p_values_combined(i))
        x = [find(unique_angles == angle_pairs(i, 1)), find(unique_angles == angle_pairs(i, 2))];
        y = [max(means) + 5, max(means) + 5];
        plot(x, y, 'k-', 'LineWidth', 1.5);
        if p_values_combined(i) < 0.05
            text(mean(x), max(means) + 7, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
        elseif p_values_combined(i) < 0.01
            text(mean(x), max(means) + 7, '**', 'FontSize', 14, 'HorizontalAlignment', 'center');
        elseif p_values_combined(i) < 0.001
            text(mean(x), max(means) + 7, '***', 'FontSize', 14, 'HorizontalAlignment', 'center');
        end
    end
end

% Add vertical lines to indicate the different angles
for i = 1:length(unique_angles)
    xline(i, '--', 'Color', colors(i,:), 'LineWidth', 1.5);
end

% Customize plot
title('Comparison of Spike Magnitudes Across Specified ZY Angle Magnitudes');
xlabel('ZY Angle Magnitude (degrees)');
ylabel('Magnitude of Spike Values (µV)');
xticks(1:length(unique_angles));
xticklabels(compose('%.1f°', unique_angles));
legend('Mean with 95% CI', 'Location', 'best');
grid on;

% Add space to the left and right of the plot
xlim([0, length(unique_angles) + 1]);

hold off;


% % Compute distances from 0
% distances_positive = abs(T.Time_Leftmost_Positive_Spike);
% distances_negative = abs(T.Time_Rightmost_Negative_Spike);

% % Correlation analysis
% [corr_pos, p_val_pos] = corr(distances_positive, T.Leftmost_Positive_Spike_Value, 'Type', 'Pearson');
% [corr_neg, p_val_neg] = corr(distances_negative, T.Rightmost_Negative_Spike_Value, 'Type', 'Pearson');

% % Display correlation results
% fprintf('Correlation Analysis:\n');
% fprintf('Positive Spikes: Correlation = %.2f, p-value = %.4f\n', corr_pos, p_val_pos);
% fprintf('Negative Spikes: Correlation = %.2f, p-value = %.4f\n', corr_neg, p_val_neg);

% % Regression analysis for positive spikes
% mdl_pos = fitlm(distances_positive, T.Leftmost_Positive_Spike_Value);
% fprintf('\nRegression Analysis for Positive Spikes:\n');
% disp(mdl_pos);

% % Regression analysis for negative spikes
% mdl_neg = fitlm(distances_negative, T.Rightmost_Negative_Spike_Value);
% fprintf('\nRegression Analysis for Negative Spikes:\n');
% disp(mdl_neg);

% % Plot the regression results
% figure;
% subplot(2,1,1);
% scatter(distances_positive, T.Leftmost_Positive_Spike_Value, 'filled');
% hold on;
% plot(mdl_pos);
% xlabel('Distance from 0 (s)');
% ylabel('Leftmost Positive Spike Value');
% title('Regression Analysis for Positive Spikes');
% legend('Data', 'Fit', 'Location', 'best');
% grid on;

% subplot(2,1,2);
% scatter(distances_negative, T.Rightmost_Negative_Spike_Value, 'filled');
% hold on;
% plot(mdl_neg);
% xlabel('Distance from 0 (s)');
% ylabel('Rightmost Negative Spike Value');
% title('Regression Analysis for Negative Spikes');
% legend('Data', 'Fit', 'Location', 'best');
% grid on;
