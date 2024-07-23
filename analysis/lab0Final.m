% Get the full path of the folder containing this script
script_folder = fileparts(mfilename('fullpath'));

% Define the path to the data folder two levels deep relative to the script folder
data_folder = fullfile(script_folder, '..', 'lab_0_data', 'converted_csvs');

% Get a list of all text files in the folder
filePattern = fullfile(data_folder, '*.txt');
txtFiles = dir(filePattern);

% Option to use the whole dataset instead of segmenting
useWholeDataset = input('Do you want to use the whole dataset instead of segmenting it? (y/n): ', 's');
useWholeDataset = strcmpi(useWholeDataset, 'y');

% Initialize arrays to store results
averageHeartRates = [];
averagePeakMagnitudes = [];
totalPeaks = [];
stdPeakHeights = [];
stdHeartRates = [];
names = {};

% Loop through each file
for k = 1:length(txtFiles)
    baseFileName = txtFiles(k).name;
    fullFileName = fullfile(data_folder, baseFileName);

    % Extract name and frequency from the filename (e.g., 'dataANDREW1000hzwristwristankle.txt')
    [tokens, ~] = regexp(baseFileName, 'data(\D+)(\d+)\D+\.txt', 'tokens', 'match');
    if ~isempty(tokens)
        name = tokens{1}{1};
        frequency = str2double(tokens{1}{2});
        
        % Check the 6th letter from the right (excluding '.txt')
        suffix_position = length(baseFileName) - 9; % -4 for '.txt' and -4 to get to the 6th from the right
        if baseFileName(suffix_position) == 't'
            name = [name '1'];
        else
            name = [name '2'];
        end
        
        names{end+1} = name;
    else
        % If frequency or name cannot be determined, skip this file
        fprintf('Skipping file %s because frequency or name is not found in the filename.\n', baseFileName);
        continue;
    end

    % Read the ECG data from the text file, skipping the header lines
    fileID = fopen(fullFileName, 'r');
    data = textscan(fileID, '%f', 'HeaderLines', 7, 'Delimiter', ',');
    fclose(fileID);
    data = data{1};

    % Check if data is empty
    if isempty(data)
        fprintf('Skipping file %s because the data is empty.\n', baseFileName);
        continue;
    end

    % Define the filter parameters
    [b, a] = butter(4, 6/(frequency/2), 'high');
    data = filter(b, a, data);
    
    % Create a time vector based on the frequency
    t = (0:length(data)-1) / frequency;
    
    % If using the whole dataset, no need to segment
    if useWholeDataset
        segmentedData = data;
        segmentedTime = t;
    else
        % Plot the ECG data
        fig1 = figure('Name', 'Full ECG Data', 'NumberTitle', 'off');
        set(fig1, 'Position', [100, 500, 600, 400]); % Position for the first figure
        plot(t, data);
        title(sprintf('%s ECG Data at %d Hz', name, frequency));
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
        segmentedData = data(idx1:idx2);
        segmentedTime = t(idx1:idx2);
    end

    % Create figure for the segmented ECG data
    fig2 = figure('Name', 'Segmented ECG Data', 'NumberTitle', 'off');
    set(fig2, 'Position', [750, 500, 600, 400]); % Position for the second figure
    plot(segmentedTime, segmentedData);
    title(sprintf('%s ECG Data at %d Hz', name, frequency));
    xlabel('Time (s)');
    ylabel('ECG Signal');
    
    % Find peaks in the segmented data
    [pks, locs] = findpeaks(segmentedData, 'MinPeakHeight', mean(segmentedData) + std(segmentedData), 'MinPeakDistance', frequency/2);
    
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
    
    title(sprintf('%s ECG Data at %d Hz with Peaks Highlighted', name, frequency));
    
    % Display results in the command window
    disp(['Name: ' name]);
    disp(['Heart Rate: ' num2str(bpm) ' BPM ± ' num2str(stdBpm)]);
    disp(['Average Peak Magnitude: ' num2str(avgPeakMagnitude) ' ± ' num2str(stdPeakHeight)]);
    disp(['Total Number of Peaks: ' num2str(length(pks))]);
end

% Create a table with the results
resultsTable0 = table(names', averageHeartRates, stdHeartRates, averagePeakMagnitudes, stdPeakHeights, totalPeaks, ...
    'VariableNames', {'Name', 'AverageHeartRate', 'StdHeartRate', 'AveragePeakMagnitude', 'StdPeakHeight', 'TotalPeaks'});

% Display summary of results
disp('Summary of Results:');
disp(resultsTable0);

% Save the table to the workspace
assignin('base', 'resultsTable0', resultsTable0);

% Unique names (excluding numbers)
uniqueNames = unique(cellfun(@(x) x(1:end-1), resultsTable0.Name, 'UniformOutput', false));

for n = 1:length(uniqueNames)
    name1 = [uniqueNames{n}, '1'];
    name2 = [uniqueNames{n}, '2'];
    
    % Extract data for each name
    data1 = resultsTable0(strcmp(resultsTable0.Name, name1), :);
    data2 = resultsTable0(strcmp(resultsTable0.Name, name2), :);
    
    % Check if data is available for both names
    if ~isempty(data1) && ~isempty(data2)
        % Calculate t-test for AverageHeartRate
        pooledStdHeartRate = sqrt(((data1.StdHeartRate.^2 .* (data1.TotalPeaks - 1)) + (data2.StdHeartRate.^2 .* (data2.TotalPeaks - 1))) / (data1.TotalPeaks + data2.TotalPeaks - 2));
        t_heartRate = (data1.AverageHeartRate - data2.AverageHeartRate) ./ (pooledStdHeartRate .* sqrt(1./data1.TotalPeaks + 1./data2.TotalPeaks));
        df_heartRate = data1.TotalPeaks + data2.TotalPeaks - 2;
        p_heartRate = 2 * tcdf(-abs(t_heartRate), df_heartRate);
        
        % Calculate t-test for AveragePeakMagnitude
        pooledStdPeakMagnitude = sqrt(((data1.StdPeakHeight.^2 .* (data1.TotalPeaks - 1)) + (data2.StdPeakHeight.^2 .* (data2.TotalPeaks - 1))) / (data1.TotalPeaks + data2.TotalPeaks - 2));
        t_peak = (data1.AveragePeakMagnitude - data2.AveragePeakMagnitude) ./ (pooledStdPeakMagnitude .* sqrt(1./data1.TotalPeaks + 1./data2.TotalPeaks));
        df_peak = data1.TotalPeaks + data2.TotalPeaks - 2;
        p_peak = 2 * tcdf(-abs(t_peak), df_peak);
        
        % Display results
        fprintf('Comparison between %s and %s:\n', name1, name2);
        fprintf('Average Heart Rate: p-value = %.4e\n', p_heartRate);
        fprintf('Average Peak Magnitude: p-value = %.4e\n\n', p_peak);
    else
        fprintf('Data not available for comparison between %s and %s\n', name1, name2);
    end
end
