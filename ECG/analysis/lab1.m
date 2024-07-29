% Get the full path of the folder containing this script
script_folder = fileparts(mfilename('fullpath'));

% Define the path to the data folder relative to the script folder
data_folder = fullfile(script_folder, '..', 'lab_1_data');

% Get a list of all text files in the data folder
data_files = dir(fullfile(data_folder, '*.txt'));

% Option to use the whole dataset instead of segmenting
useWholeDataset = input('Do you want to use the whole dataset instead of segmenting it? (y/n): ', 's');
useWholeDataset = strcmpi(useWholeDataset, 'y');

% Define the sampling rate (samples per second)
sampling_rate = 1000; % Adjust this value based on your specific needs

% Initialize arrays to store results
averageHeartRates = [];
averagePeakMagnitudes = [];
totalPeaks = [];
stdPeakHeights = [];
stdHeartRates = [];
names = {};

% Define the filter parameters
[b, a] = butter(4, 6/(sampling_rate/2), 'high');

% Loop through each file and process it
for k = 1:length(data_files)
    baseFileName = data_files(k).name;
    fullFileName = fullfile(data_folder, baseFileName);

    % Extract name from the filename
    name = regexp(baseFileName, 'Lab1(\w+)', 'tokens', 'once');
    if ~isempty(name)
        name = name{1};
        names{end+1} = name;
    else
        name = 'Unknown';
        names{end+1} = name;
    end

    % Read the data from the text file
    data = readmatrix(fullFileName);
    
    % Check if the file is inverted and flip the data if necessary
    if contains(baseFileName, 'Inverted')
        data = -data;
    end
    
    % Create a time vector based on the sampling rate
    t = (0:length(data)-1)' / sampling_rate;
    
    % Apply the high-pass filter to the data
    filtered_data = filter(b, a, data);
    
    % If using the whole dataset, no need to segment
    if useWholeDataset
        segmentedData = filtered_data;
        segmentedTime = t;
    else
        % Plot the ECG data
        fig1 = figure('Name', 'Full ECG Data', 'NumberTitle', 'off');
        set(fig1, 'Position', [100, 500, 600, 400]);
        plot(t, filtered_data);
        title(sprintf('%s ECG Data', name));
        xlabel('Time (s)');
        ylabel('ECG Signal');
        
        % Capture two clicks from the user
        disp('Please click two points on the graph to select the segment.');
        [x, ~] = ginput(2);
        x1 = x(1);
        x2 = x(2);
        
        % Close the original plot after segmenting
        close(fig1);
        
        % Find the indices corresponding to the clicked x values
        idx1 = find(t >= x1, 1, 'first');
        idx2 = find(t <= x2, 1, 'last');
        
        % Segment the data
        segmentedData = filtered_data(idx1:idx2);
        segmentedTime = t(idx1:idx2);
    end
    
    % Create figure for the segmented ECG data
    fig2 = figure('Name', 'Segmented ECG Data', 'NumberTitle', 'off');
    set(fig2, 'Position', [750, 500, 600, 400]);
    plot(segmentedTime, segmentedData);
    title(sprintf('%s ECG Data', name));
    xlabel('Time (s)');
    ylabel('ECG Signal');
    
    % Find peaks in the segmented data
    [pks, locs] = findpeaks(segmentedData, 'MinPeakHeight', mean(segmentedData) + std(segmentedData), 'MinPeakDistance', sampling_rate/2);
    
    % Check if peaks were found
    if isempty(pks)
        fprintf('No peaks found in file %s.\n', baseFileName);
        continue;
    end
    
    % Calculate the time between peaks
    peakIntervals = diff(segmentedTime(locs));
    
    % Calculate BPM and its standard deviation
    avgPeakInterval = mean(peakIntervals);
    bpm = 60 / avgPeakInterval;
    stdBpm = std(60 ./ peakIntervals);
    
    % Calculate average peak magnitude and its standard deviation
    avgPeakMagnitude = mean(pks);
    stdPeakHeight = std(pks);
    
    % Store results
    averageHeartRates = [averageHeartRates; bpm];
    averagePeakMagnitudes = [averagePeakMagnitudes; avgPeakMagnitude];
    totalPeaks = [totalPeaks; length(pks)];
    stdPeakHeights = [stdPeakHeights; stdPeakHeight];
    stdHeartRates = [stdHeartRates; stdBpm];
    
    % Highlight the peaks on the plot
    hold on;
    plot(segmentedTime(locs), segmentedData(locs), 'ro');
    hold off;
    
    % Display BPM, average peak magnitude, standard deviations, and total peaks on the graph
    text(mean(segmentedTime), max(segmentedData), sprintf('BPM: %.2f ± %.2f\nAvg Peak Magnitude: %.2f ± %.2f\nTotal Peaks: %d', ...
        bpm, stdBpm, avgPeakMagnitude, stdPeakHeight, length(pks)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'BackgroundColor', 'w');
    
    title(sprintf('%s ECG Data with Peaks Highlighted', name));
    
    % Display results in the command window
    disp(['Name: ' name]);
    disp(['Heart Rate: ' num2str(bpm) ' BPM ± ' num2str(stdBpm)]);
    disp(['Average Peak Magnitude: ' num2str(avgPeakMagnitude) ' ± ' num2str(stdPeakHeight)]);
    disp(['Total Number of Peaks: ' num2str(length(pks))]);
end

% Create a table with the results
resultsTable = table(names', averageHeartRates, stdHeartRates, averagePeakMagnitudes, stdPeakHeights, totalPeaks, ...
    'VariableNames', {'Name', 'AverageHeartRate', 'StdHeartRate', 'AveragePeakMagnitude', 'StdPeakHeight', 'TotalPeaks'});

% Display summary of results
disp('Summary of Results:');
disp(resultsTable);

% Save the table to the workspace
assignin('base', 'resultsTable', resultsTable);

% Extract data for the inverted and normal datasets
invertedData = resultsTable(contains(resultsTable.Name, 'Inverted'), :);
normalData = resultsTable(contains(resultsTable.Name, 'Normal'), :);

% Perform t-tests between inverted and normal datasets
n1 = invertedData.TotalPeaks;
n2 = normalData.TotalPeaks;

% Calculate pooled standard deviation for heart rate
sp_heart_rate = sqrt(((n1 - 1) .* invertedData.StdHeartRate.^2 + (n2 - 1) .* normalData.StdHeartRate.^2) / (n1 + n2 - 2));
t_heart_rate = (invertedData.AverageHeartRate - normalData.AverageHeartRate) ./ (sp_heart_rate .* sqrt(1/n1 + 1/n2));
df_heart_rate = n1 + n2 - 2;
p_heart_rate = 2 * tcdf(-abs(t_heart_rate), df_heart_rate);
h_heart_rate = p_heart_rate < 0.05;

% Calculate pooled standard deviation for peak size
sp_peak_size = sqrt(((n1 - 1) .* invertedData.StdPeakHeight.^2 + (n2 - 1) .* normalData.StdPeakHeight.^2) / (n1 + n2 - 2));
t_peak_size = (invertedData.AveragePeakMagnitude - normalData.AveragePeakMagnitude) ./ (sp_peak_size .* sqrt(1/n1 + 1/n2));
df_peak_size = n1 + n2 - 2;
p_peak_size = 2 * tcdf(-abs(t_peak_size), df_peak_size);
h_peak_size = p_peak_size < 0.05;

% Display t-test results
disp('T-test Results between Inverted and Normal datasets:');
fprintf('Heart Rate: h = %d, p = %.4f\n', h_heart_rate, p_heart_rate);
fprintf('Peak Size: h = %d, p = %.4f\n', h_peak_size, p_peak_size);

% Save the t-test results to the workspace
assignin('base', 'tTestResults_HeartRate', table(h_heart_rate, p_heart_rate, 'VariableNames', {'h', 'p'}));
assignin('base', 'tTestResults_PeakSize', table(h_peak_size, p_peak_size, 'VariableNames', {'h', 'p'}));
