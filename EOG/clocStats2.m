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

% Loop through each file and process the data
for k = 1:length(file_names)
    file_name = file_names{k};
    data = data_map(file_name);

    % Process the data using the process_data function
    % Assuming process_data returns processed data
    calibrated_data = process_data(data, 200, 2, file_name);

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

    % Remove specified peaks for RightStartCounterClockwise.txt
    if strcmp(file_name, 'RightStartCounterClockwise.txt')
        remove_indices = [2, 4, 9];
        if length(significant_locs_positive) >= max(remove_indices)
            significant_locs_positive(remove_indices) = [];
        end
    end

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
end

% Display the results
disp(results);

% Perform linear regression between specified pairs of points
regression_results = table();

% Identify indices for pairs of points for regression
positive_indices = find(strcmp(results.Peak_Type, 'Positive'));
negative_indices = find(strcmp(results.Peak_Type, 'Negative'));

% Ensure there are enough peaks to perform the requested regressions
num_regressions = min(floor(length(positive_indices)), floor(length(negative_indices)/2));

for i = 1:num_regressions
    % First linear regression
    neg_idx1 = negative_indices(2*i-1);
    pos_idx1 = positive_indices(i);
    X1 = [results.Time(neg_idx1), results.Time(pos_idx1)];
    Y1 = [results.Peak_Value(neg_idx1), results.Peak_Value(pos_idx1)];
    
    % Fit linear model
    p1 = polyfit(X1, Y1, 1);
    yfit1 = polyval(p1, X1);
    
    % Calculate R-squared value
    yresid1 = Y1 - yfit1;
    SSresid1 = sum(yresid1.^2);
    SStotal1 = (length(Y1)-1) * var(Y1);
    rsq1 = 1 - SSresid1/SStotal1;

    % Store regression results
    T_regression1 = table(X1', Y1', yfit1', repmat({results.Filename{neg_idx1}}, 2, 1), repmat(p1(1), 2, 1), repmat(p1(2), 2, 1), repmat(rsq1, 2, 1), ...
                         'VariableNames', {'Time', 'Peak_Value', 'Fit', 'Filename', 'Slope', 'Intercept', 'R_squared'});
    regression_results = [regression_results; T_regression1];
    
    % Second linear regression
    pos_idx2 = positive_indices(i);
    neg_idx2 = negative_indices(2*i);
    X2 = [results.Time(pos_idx2), results.Time(neg_idx2)];
    Y2 = [results.Peak_Value(pos_idx2), results.Peak_Value(neg_idx2)];
    
    % Fit linear model
    p2 = polyfit(X2, Y2, 1);
    yfit2 = polyval(p2, X2);
    
    % Calculate R-squared value
    yresid2 = Y2 - yfit2;
    SSresid2 = sum(yresid2.^2);
    SStotal2 = (length(Y2)-1) * var(Y2);
    rsq2 = 1 - SSresid2/SStotal2;

    % Store regression results
    T_regression2 = table(X2', Y2', yfit2', repmat({results.Filename{pos_idx2}}, 2, 1), repmat(p2(1), 2, 1), repmat(p2(2), 2, 1), repmat(rsq2, 2, 1), ...
                         'VariableNames', {'Time', 'Peak_Value', 'Fit', 'Filename', 'Slope', 'Intercept', 'R_squared'});
    regression_results = [regression_results; T_regression2];
end

% Display the regression results
disp(regression_results);

% Plot the data and regression lines for visualization
figure;
num_subplots = num_regressions * 2;
subplot_idx = 1;

for i = 1:num_regressions
    % First linear regression
    neg_idx1 = negative_indices(2*i-1);
    pos_idx1 = positive_indices(i);
    X1 = [results.Time(neg_idx1), results.Time(pos_idx1)];
    Y1 = [results.Peak_Value(neg_idx1), results.Peak_Value(pos_idx1)];
    
    % Fit linear model
    p1 = polyfit(X1, Y1, 1);
    yfit1 = polyval(p1, X1);
    
    % Calculate R-squared value
    yresid1 = Y1 - yfit1;
    SSresid1 = sum(yresid1.^2);
    SStotal1 = (length(Y1)-1) * var(Y1);
    rsq1 = 1 - SSresid1/SStotal1;

    % Plot the first linear regression
    subplot(num_subplots, 1, subplot_idx);
    plot(X1, Y1, 'o', 'DisplayName', 'Data Points');
    hold on;
    plot(X1, yfit1, '-', 'LineWidth', 1.5, 'DisplayName', 'Fit Line');
    title(sprintf('Regression Line %d: Slope: %.2f, Intercept: %.2f, R^2: %.2f', subplot_idx, p1(1), p1(2), rsq1));
    xlabel('Time (s)');
    ylabel('EOG Signal');
    legend('show');
    hold off;
    
    subplot_idx = subplot_idx + 1;
    
    % Second linear regression
    pos_idx2 = positive_indices(i);
    neg_idx2 = negative_indices(2*i);
    X2 = [results.Time(pos_idx2), results.Time(neg_idx2)];
    Y2 = [results.Peak_Value(pos_idx2), results.Peak_Value(neg_idx2)];
    
    % Fit linear model
    p2 = polyfit(X2, Y2, 1);
    yfit2 = polyval(p2, X2);
    
    % Calculate R-squared value
    yresid2 = Y2 - yfit2;
    SSresid2 = sum(yresid2.^2);
    SStotal2 = (length(Y2)-1) * var(Y2);
    rsq2 = 1 - SSresid2/SStotal2;

    % Plot the second linear regression
    subplot(num_subplots, 1, subplot_idx);
    plot(X2, Y2, 'o', 'DisplayName', 'Data Points');
    hold on;
    plot(X2, yfit2, '-', 'LineWidth', 1.5, 'DisplayName', 'Fit Line');
    title(sprintf('Regression Line %d: Slope: %.2f, Intercept: %.2f, R^2: %.2f', subplot_idx, p2(1), p2(2), rsq2));
    xlabel('Time (s)');
    ylabel('EOG Signal');
    legend('show');
    hold off;
    
    subplot_idx = subplot_idx + 1;
end

sgtitle('Individual Regression Lines with Statistics');



