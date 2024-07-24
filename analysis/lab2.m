% Get the full path of the folder containing this script
script_folder = fileparts(mfilename('fullpath'));

% Define the path to the data folder relative to the script folder
data_folder = fullfile(script_folder, '..', 'lab_2_data');

% Get a list of all text files in the folder
filePattern = fullfile(data_folder, '*.txt');
txtFiles = dir(filePattern);

% Option to use the whole dataset instead of segmenting
useWholeDataset = input('Do you want to use the whole dataset instead of segmenting it? (y/n): ', 's');
useWholeDataset = strcmpi(useWholeDataset, 'y');

% Initialize arrays to store results
averagePeakHeights = [];
stdPeakHeights = [];
frequencies = [];
bpms = [];
stdBPMs = [];  % Change the name to stdBPMs
totalPeaks = [];
names = {};

% Loop through each file
for k = 1:length(txtFiles)
    baseFileName = txtFiles(k).name;
    fullFileName = fullfile(data_folder, baseFileName);

    % Extract frequency from the filename (e.g., 'Andrew1000hzWAA.txt')
    freq = regexp(baseFileName, '\d+', 'match');
    if ~isempty(freq)
        frequency = str2double(freq{1});
        frequencies = [frequencies; frequency];
    else
        % If frequency cannot be determined, skip this file
        fprintf('Skipping file %s because frequency is not found in the filename.\n', baseFileName);
        continue;
    end

    % Extract name from the filename
    name = regexp(baseFileName, 'data(\w+)\d+hz', 'tokens', 'once');
    if ~isempty(name)
        name = name{1};
        names{end+1} = name;
    else
        name = 'Unknown';
        names{end+1} = name;
    end

    % Define the filter parameters
    [b, a] = butter(4, 6/(frequency/2), 'high');
    
    % Read the ECG data from the text file
    fileID = fopen(fullFileName, 'r');
    data = fscanf(fileID, '%f,');
    fclose(fileID);

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
        title(sprintf('ECG Data at %d Hz', frequency));
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
    title(sprintf('Segmented ECG Data at %d Hz', frequency));
    xlabel('Time (s)');
    ylabel('ECG Signal');
    
    % Find peaks in the segmented data
    [pks, locs] = findpeaks(segmentedData, 'MinPeakHeight', mean(segmentedData) + std(segmentedData), 'MinPeakDistance', frequency/2);
    
    % Calculate the time between peaks
    peakIntervals = diff(segmentedTime(locs));
    
    % Calculate BPM for each interval
    bpms_temp = 60 ./ peakIntervals;
    bpm = mean(bpms_temp);
    stdBPM = std(bpms_temp);  % Calculate the standard deviation of BPMs
    
    % Calculate average and standard deviation of peak heights
    avgPeakHeight = mean(pks);
    stdPeakHeight = std(pks);
    
    % Store results
    averagePeakHeights = [averagePeakHeights; avgPeakHeight];
    stdPeakHeights = [stdPeakHeights; stdPeakHeight];
    bpms = [bpms; bpm];
    stdBPMs = [stdBPMs; stdBPM];  % Store the standard deviation of BPMs
    totalPeaks = [totalPeaks; length(pks)];
    
    % Highlight the peaks on the plot
    hold on;
    plot(segmentedTime(locs), segmentedData(locs), 'ro');
    hold off;
    
    % Display BPM, standard deviation, average peak height, and total peaks on the graph
    text(mean(segmentedTime), max(segmentedData), sprintf('BPM: %.2f\nStd Dev BPM: %.2f\nAvg Peak Height: %.2f\nStd Dev Height: %.2f\nTotal Peaks: %d', ...
        bpm, stdBPM, avgPeakHeight, stdPeakHeight, length(pks)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'BackgroundColor', 'w');
    
    title(sprintf('Segmented ECG Data at %d Hz with Peaks Highlighted', frequency));
    
    % Display BPM and standard deviation in the command window
    disp(['Heart Rate: ' num2str(bpm) ' BPM']);
    disp(['Standard Deviation of BPMs: ' num2str(stdBPM)]);
    disp(['Average Peak Height: ' num2str(avgPeakHeight)]);
    disp(['Standard Deviation of Peak Heights: ' num2str(stdPeakHeight)]);
    disp(['Total Number of Peaks: ' num2str(length(pks))]);
end

% Normalize the data
normalizedAvgPeakHeights = (averagePeakHeights - min(averagePeakHeights)) / (max(averagePeakHeights) - min(averagePeakHeights));
normalizedStdPeakHeights = (stdPeakHeights - min(stdPeakHeights)) / (max(stdPeakHeights) - min(stdPeakHeights));

% Flip the normalized x-axis values for Euclidean distance calculation
flippedNormalizedAvgPeakHeights = 1 - normalizedAvgPeakHeights;

% Calculate the Euclidean distance from the origin (0,0) for each point with the flipped x-axis
euclideanDistances = sqrt(flippedNormalizedAvgPeakHeights.^2 + normalizedStdPeakHeights.^2);

% Plot normalized average peak height against normalized standard deviation of peak heights
figure('Name', 'Normalized Average Peak Height vs Normalized Std Dev of Peak Heights', 'NumberTitle', 'off');
scatter(normalizedAvgPeakHeights, normalizedStdPeakHeights);
xlabel('Normalized Average Peak Height');
ylabel('Normalized Standard Deviation of Peak Heights');
title('Normalized Average Peak Height vs Normalized Standard Deviation of Peak Heights');
set(gca, 'XDir', 'reverse'); % Flip the x-axis

% Label the points with the frequency (Hz) and Euclidean distance
for i = 1:length(frequencies)
    text(normalizedAvgPeakHeights(i), normalizedStdPeakHeights(i), ...
        sprintf('%d Hz\nDist: %.2f', frequencies(i), euclideanDistances(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Create a table with the results
names(:) = {'ANDREW2'}; % Set all names to 'ANDREW2'
resultsTable2 = table(names', frequencies, bpms, stdBPMs, averagePeakHeights, stdPeakHeights, totalPeaks, ...
    'VariableNames', {'Name', 'Frequency', 'BPM', 'StdBPM', 'AveragePeakHeight', 'StdPeakHeight', 'TotalPeaks'});

% Display summary of results
disp('Summary of Results:');
disp(resultsTable2);

% Save the table to the workspace
assignin('base', 'resultsTable2', resultsTable2);

% Sort the results table by frequency in descending order
sortedResultsTable = sortrows(resultsTable2, 'Frequency', 'descend');

% Initialize arrays to store results
frequencies1 = [];
frequencies2 = [];
pValuesBPM = [];
pValuesPeakHeight = [];
hValuesBPM = [];
hValuesPeakHeight = [];

% Loop through each row in the sorted results table and compare with the next lower frequency dataset
for i = 1:height(sortedResultsTable) - 1
    % Extract data for the current dataset
    M1_BPM = sortedResultsTable.BPM(i);
    S1_BPM = sortedResultsTable.StdBPM(i);
    n1 = sortedResultsTable.TotalPeaks(i);

    M1_PeakHeight = sortedResultsTable.AveragePeakHeight(i);
    S1_PeakHeight = sortedResultsTable.StdPeakHeight(i);

    % Extract data for the next lower frequency dataset
    M2_BPM = sortedResultsTable.BPM(i + 1);
    S2_BPM = sortedResultsTable.StdBPM(i + 1);
    n2 = sortedResultsTable.TotalPeaks(i + 1);

    M2_PeakHeight = sortedResultsTable.AveragePeakHeight(i + 1);
    S2_PeakHeight = sortedResultsTable.StdPeakHeight(i + 1);

    % Calculate t-statistic for BPM
    t_BPM = (M1_BPM - M2_BPM) / sqrt((S1_BPM^2 / n1) + (S2_BPM^2 / n2));
    % Calculate degrees of freedom for BPM
    df_BPM = ((S1_BPM^2 / n1) + (S2_BPM^2 / n2))^2 / ((S1_BPM^2 / n1)^2 / (n1 - 1) + (S2_BPM^2 / n2)^2 / (n2 - 1));
    % Calculate p-value for BPM (two-tailed test)
    p_BPM = max(2 * (1 - tcdf(abs(t_BPM), df_BPM)), eps); % Ensure p-value is not exactly 0
    % Calculate h value for BPM (reject null hypothesis if p-value < 0.05)
    h_BPM = p_BPM < 0.05;

    % Calculate t-statistic for Peak Height (one-tailed test)
    t_PeakHeight = (M1_PeakHeight - M2_PeakHeight) / sqrt((S1_PeakHeight^2 / n1) + (S2_PeakHeight^2 / n2));
    % Calculate degrees of freedom for Peak Height
    df_PeakHeight = ((S1_PeakHeight^2 / n1) + (S2_PeakHeight^2 / n2))^2 / ((S1_PeakHeight^2 / n1)^2 / (n1 - 1) + (S2_PeakHeight^2 / n2)^2 / (n2 - 1));
    % Calculate p-value for Peak Height (one-tailed test)
    p_PeakHeight = max(1 - tcdf(t_PeakHeight, df_PeakHeight), eps); % Ensure p-value is not exactly 0
    % Calculate h value for Peak Height (reject null hypothesis if p-value < 0.05)
    h_PeakHeight = p_PeakHeight < 0.05;

    % Store results
    frequencies1 = [frequencies1; sortedResultsTable.Frequency(i)];
    frequencies2 = [frequencies2; sortedResultsTable.Frequency(i + 1)];
    pValuesBPM = [pValuesBPM; p_BPM];
    pValuesPeakHeight = [pValuesPeakHeight; p_PeakHeight];
    hValuesBPM = [hValuesBPM; h_BPM];
    hValuesPeakHeight = [hValuesPeakHeight; h_PeakHeight];

    % Display results
    fprintf('Comparison of %d Hz dataset with %d Hz dataset:\n', sortedResultsTable.Frequency(i), sortedResultsTable.Frequency(i + 1));
    fprintf('  P-value for BPM: %.4e (h = %d)\n', p_BPM, h_BPM);
    fprintf('  P-value for Average Peak Height: %.4e (h = %d)\n\n', p_PeakHeight, h_PeakHeight);
end

% Create summary table of p-values and h-values
comparisonTable = table(frequencies1, frequencies2, pValuesBPM, hValuesBPM, pValuesPeakHeight, hValuesPeakHeight, ...
    'VariableNames', {'Frequency1', 'Frequency2', 'PValueBPM', 'HBPM', 'PValuePeakHeight', 'HPeakHeight'});

% Extract and sort data for the scatter plot
[sortedFrequencies, sortIdx] = sort(resultsTable2.Frequency, 'descend');
sortedAveragePeakHeights = resultsTable2.AveragePeakHeight(sortIdx);
sortedStdPeakHeights = resultsTable2.StdPeakHeight(sortIdx);

% Create a categorical array for frequencies, ensuring they are in the desired order
frequencies = categorical(string(sortedFrequencies), string(sortedFrequencies), 'Ordinal', true);

% Create the scatter plot with error bars
figure;
errorbar(frequencies, sortedAveragePeakHeights, sortedStdPeakHeights, 's', 'MarkerSize', 10, 'LineWidth', 1.5);
title('Scatter Plot of Peak Height with Error Bars');
xlabel('Frequency (Hz)');
ylabel('Average Peak Height');
grid on;

% Annotate the plot with relevant information
for i = 1:length(frequencies)
    text(double(frequencies(i)) - 0.1, sortedAveragePeakHeights(i), sprintf('%.2f ± %.2f', sortedAveragePeakHeights(i), sortedStdPeakHeights(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
end
