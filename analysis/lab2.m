% Get the full path of the folder containing this script
%script_folder = fileparts(mfilename('fullpath'));

% Define the path to the data folder relative to the script folder
data_folder = fullfile(script_folder, '../lab_2_data/');

% Get a list of all text files in the folder
filePattern = fullfile(data_folder, '*.txt');
txtFiles = dir(filePattern);

% Option to use the whole dataset instead of segmenting
useWholeDataset = input('Do you want to use the whole dataset instead of segmenting it? (y/n): ', 's');
useWholeDataset = strcmpi(useWholeDataset, 'y');

% Initialize arrays to store average peak height and standard deviation
averagePeakHeights = [];
stdPeakHeights = [];
frequencies = [];

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
    
    % Calculate BPM and standard deviation
    avgPeakInterval = mean(peakIntervals);
    bpm = 60 / avgPeakInterval;
    stdPeakInterval = std(peakIntervals);
    
    % Calculate average and standard deviation of peak heights
    avgPeakHeight = mean(pks);
    stdPeakHeight = std(pks);
    
    % Store average peak height and standard deviation
    averagePeakHeights = [averagePeakHeights; avgPeakHeight];
    stdPeakHeights = [stdPeakHeights; stdPeakHeight];
    
    % Highlight the peaks on the plot
    hold on;
    plot(segmentedTime(locs), segmentedData(locs), 'ro');
    hold off;
    
    % Display BPM, standard deviation, average peak height, and total peaks on the graph
    text(mean(segmentedTime), max(segmentedData), sprintf('BPM: %.2f\nStd Dev Interval: %.2f\nAvg Peak Height: %.2f\nStd Dev Height: %.2f\nTotal Peaks: %d', ...
        bpm, stdPeakInterval, avgPeakHeight, stdPeakHeight, length(pks)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'BackgroundColor', 'w');
    
    title(sprintf('Segmented ECG Data at %d Hz with Peaks Highlighted', frequency));
    
    % Display BPM and standard deviation in the command window
    disp(['Heart Rate: ' num2str(bpm) ' BPM']);
    disp(['Standard Deviation of Peak Intervals: ' num2str(stdPeakInterval)]);
    disp(['Average Peak Height: ' num2str(avgPeakHeight)]);
    disp(['Standard Deviation of Peak Heights: ' num2str(stdPeakHeight)]);
    disp(['Total Number of Peaks: ' num2str(length(pks))]);
    
end

% Plot average peak height against standard deviation of peak heights
figure('Name', 'Average Peak Height vs Std Dev of Peak Heights', 'NumberTitle', 'off');
scatter(averagePeakHeights, stdPeakHeights);
xlabel('Average Peak Height');
ylabel('Standard Deviation of Peak Heights');
title('Average Peak Height vs Standard Deviation of Peak Heights');

% Label the points with the frequency (Hz) of the data point
for i = 1:length(frequencies)
    text(averagePeakHeights(i), stdPeakHeights(i), sprintf('%d Hz', frequencies(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
