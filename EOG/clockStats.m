% Load the table containing the calibration data
load('calibration2_data_table.mat');

% Extract data for the specified files
file_names = {'RightStartCounterClockwise.txt', 'rightStartClockwise.txt'};
data_map = containers.Map;

for i = 1:height(data_table)
    if ismember(data_table.Filename{i}, file_names)
        data_map(data_table.Filename{i}) = data_table.Data{i};
    end
end

% Check if data is found
if isempty(data_map)
    error('Specified data files not found in data_table.');
end

% Initialize arrays to store the results
results = table();
linear_regression_results = {};
segment_counter = 0;

% Initialize figures for subplots and overlay
figure1 = figure;
figure2 = figure;

subplot_idx_1 = 1;  % Subplot index for figure 1
subplot_idx_2 = 1;  % Subplot index for figure 2

% Loop through each file and process the data
for k = 1:length(file_names)
    file_name = file_names{k};
    data = data_map(file_name);

    % Process the data using the process_data function
    % Assuming process_data returns processed data
    % Here you need to implement or call your `process_data` function correctly
    % calibrated_data = process_data(data, 200, 2, file_name);
    % For demonstration purposes, let's assume the processed data is the same as the input data
    calibrated_data = data;

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
    significant_locs_negative = locs_negative(pks_negative > -threshold_negative);

    % Debugging statements
    disp(['Processing ', file_name]);
    disp(['Positive peaks found: ', num2str(length(significant_locs_positive))]);
    disp(['Negative peaks found: ', num2str(length(significant_locs_negative))]);

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
    else
        grouped_locs_positive = {};
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
    else
        grouped_locs_negative = {};
    end

    % Get the highest positive spikes
    highest_positive_spikes = [];
    for i = 1:length(grouped_locs_positive)
        group = grouped_locs_positive{i};
        [~, max_idx] = max(calibrated_data(group));
        highest_positive_spikes = [highest_positive_spikes; group(max_idx)];
    end

    % Get the lowest negative spikes
    lowest_negative_spikes = [];
    for i = 1:length(grouped_locs_negative)
        group = grouped_locs_negative{i};
        [~, min_idx] = min(calibrated_data(group));
        lowest_negative_spikes = [lowest_negative_spikes; group(min_idx)];
    end

    % Remove specified peaks for RightStartCounterClockwise.txt
    if strcmp(file_name, 'RightStartCounterClockwise.txt')
        remove_indices_positive = [2];
        remove_indices_negative = [4, 9];
        highest_positive_spikes(remove_indices_positive) = [];
        lowest_negative_spikes(remove_indices_negative) = [];
    end

    % Convert indices to time
    time_highest_positive_spikes = highest_positive_spikes / 200;
    time_lowest_negative_spikes = lowest_negative_spikes / 200;

    % Create tables for highest and lowest peaks
    T_positive = table(time_highest_positive_spikes, calibrated_data(highest_positive_spikes), repmat({file_name}, length(time_highest_positive_spikes), 1), ...
              'VariableNames', {'Time', 'Peak_Value', 'Filename'});
    T_negative = table(time_lowest_negative_spikes, calibrated_data(lowest_negative_spikes), repmat({file_name}, length(time_lowest_negative_spikes), 1), ...
              'VariableNames', {'Time', 'Peak_Value', 'Filename'});
    
    % Add an indicator for peak type
    T_positive.Peak_Type = repmat({'Positive'}, length(time_highest_positive_spikes), 1);
    T_negative.Peak_Type = repmat({'Negative'}, length(time_lowest_negative_spikes), 1);
    
    % Concatenate the results
    results = [results; T_positive; T_negative];

    % Perform linear regressions and plot the regression lines
    segments = [
        1, 1;  % Segment 1: Negative 1 to Positive 1
        2, 1;  % Segment 2: Negative 2 to Positive 1
        3, 2;  % Segment 3: Negative 3 to Positive 2
        4, 2;  % Segment 4: Negative 4 to Positive 2
        5, 3;  % Segment 5: Negative 5 to Positive 3
        6, 3;  % Segment 6: Negative 6 to Positive 3
        7, 4;  % Segment 7: Negative 7 to Positive 4
        8, 4;  % Segment 8: Negative 8 to Positive 4
    ];
    
    for j = 1:size(segments, 1)
        neg_idx = segments(j, 1);
        pos_idx = segments(j, 2);

        if neg_idx <= length(lowest_negative_spikes) && pos_idx <= length(highest_positive_spikes)
            start_idx = lowest_negative_spikes(neg_idx);
            end_idx = highest_positive_spikes(pos_idx);
            segment_counter = segment_counter + 1;
            
            if mod(segment_counter, 2) == 1
                segment_type = sprintf('0° - 180° (%d)', (segment_counter + 1) / 2);
            else
                segment_type = sprintf('180° - 360° (%d)', segment_counter / 2);
            end
            
            if start_idx < end_idx
                x1 = (start_idx:end_idx)' / 200;
                y1 = calibrated_data(start_idx:end_idx);
            else
                x1 = (end_idx:start_idx)' / 200;
                y1 = calibrated_data(end_idx:start_idx);
            end

            % Perform linear regression using fitlm
            lm = fitlm(x1, y1);
            y1_fit = predict(lm, x1);
            
            % Confidence intervals for the fitted values
            [y1_fit, y1_ci] = predict(lm, x1, 'Alpha', 0.05, 'Prediction', 'observation');
            
            % Store regression results
            p1 = lm.Coefficients.Estimate';
            r_squared = lm.Rsquared.Ordinary;
            p_value = lm.Coefficients.pValue(2); % p-value for the slope
            confidence_intervals = lm.coefCI;
            
            % Perform t-test on the slope
            slope = lm.Coefficients.Estimate(2);
            t_stat = (slope - 0) / lm.Coefficients.SE(2);
            df = lm.DFE;
            t_critical = tinv(0.975, df);
            p_val = 2 * (1 - tcdf(abs(t_stat), df));
            
            % Add results of t-test to the regression table
            linear_regression_results = [linear_regression_results; {file_name, segment_type, slope, r_squared, p_value, confidence_intervals, t_stat, p_val}];
            
            % Create subplot in figure 1 for RightStartCounterClockwise.txt
            if strcmp(file_name, 'RightStartCounterClockwise.txt')
                figure(figure1);
                subplot(2, 4, subplot_idx_1);
                plot(x1, y1, 'bo'); % Original data points used in regression
                hold on;
                plot(x1, y1_fit, 'r-'); % Regression line
                fill([x1; flipud(x1)], [y1_ci(:,1); flipud(y1_ci(:,2))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % Confidence intervals
                title(segment_type);
                xlabel('Time (s)');
                ylabel('Data Value');
                legend('Data Points', 'Fit', 'Confidence Interval');
                text(0.05, 0.95, sprintf('R^2 = %.2f\np = %.2e\nCI = [%.2f, %.2f]\nt = %.2f\np-val = %.2e', r_squared, p_value, confidence_intervals(2, 1), confidence_intervals(2, 2), t_stat, p_val), 'Units', 'normalized', 'VerticalAlignment', 'top');
                hold off;
                subplot_idx_1 = subplot_idx_1 + 1;
            % Create subplot in figure 2 for rightStartClockwise.txt
            else
                figure(figure2);
                subplot(2, 4, subplot_idx_2);
                plot(x1, y1, 'bo'); % Original data points used in regression
                hold on;
                plot(x1, y1_fit, 'r-'); % Regression line
                fill([x1; flipud(x1)], [y1_ci(:,1); flipud(y1_ci(:,2))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % Confidence intervals
                title(segment_type);
                xlabel('Time (s)');
                ylabel('Data Value');
                legend('Data Points', 'Fit', 'Confidence Interval');
                text(0.05, 0.95, sprintf('R^2 = %.2f\np = %.2e\nCI = [%.2f, %.2f]\nt = %.2f\np-val = %.2e', r_squared, p_value, confidence_intervals(2, 1), confidence_intervals(2, 2), t_stat, p_val), 'Units', 'normalized', 'VerticalAlignment', 'top');
                hold off;
                subplot_idx_2 = subplot_idx_2 + 1;
            end
        end
    end
    
    % Plot original data with all linear regressions overlayed
    figure;
    hold on;
    plot((1:numel(calibrated_data))/200, calibrated_data, 'k'); % Full data
    
    for j = 1:size(segments, 1)
        neg_idx = segments(j, 1);
        pos_idx = segments(j, 2);

        if neg_idx <= length(lowest_negative_spikes) && pos_idx <= length(highest_positive_spikes)
            start_idx = lowest_negative_spikes(neg_idx);
            end_idx = highest_positive_spikes(pos_idx);
            segment_counter = segment_counter + 1;
            
            if mod(segment_counter, 2) == 1
                segment_type = sprintf('0° - 180° (%d)', (segment_counter + 1) / 2);
            else
                segment_type = sprintf('180° - 360° (%d)', segment_counter / 2);
            end
            
            if start_idx < end_idx
                x1 = (start_idx:end_idx)' / 200;
                y1 = calibrated_data(start_idx:end_idx);
            else
                x1 = (end_idx:start_idx)' / 200;
                y1 = calibrated_data(end_idx:start_idx);
            end

            % Perform linear regression using fitlm
            lm = fitlm(x1, y1);
            y1_fit = predict(lm, x1);
            
            % Confidence intervals for the fitted values
            [y1_fit, y1_ci] = predict(lm, x1, 'Alpha', 0.05, 'Prediction', 'observation');
            
            % Overlay regression line on the original data
            plot(x1, y1_fit, 'r-');
            fill([x1; flipud(x1)], [y1_ci(:,1); flipud(y1_ci(:,2))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % Confidence intervals
        end
    end
    
    hold off;
    title([file_name, ' - Data with Linear Regressions']);
    xlabel('Time (s)');
    ylabel('Data Value');
    legend('Data', 'Linear Regressions', 'Confidence Interval');
end

% Display the results
disp(results);

% Display linear regression results
linear_regression_table = cell2table(linear_regression_results, 'VariableNames', {'Filename', 'Segment_Type', 'Slope', 'R_squared', 'p_value', 'Confidence_Intervals', 't_stat', 'p_val'});
disp(linear_regression_table);
