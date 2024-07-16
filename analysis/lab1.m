% Get the full path of the folder containing this script
script_folder = fileparts(mfilename('fullpath'));

% Define the path to the data folder relative to the script folder
data_folder = fullfile(script_folder, '../lab_1_data/');

% Get a list of all CSV files in the data folder
data_files = dir(fullfile(data_folder, '*.txt'));

% Define the sampling rate (samples per second)
% Adjust this value based on your specific needs
sampling_rate = 1000; % example: 1000 samples per second

% Initialize an empty table
combined_data = table;

% Define the filter parameters
[b, a] = butter(4, 6/(sampling_rate/2), 'high');

% Loop through each file and process it
for k = 1:length(data_files)
    % Get the full path of the current file
    filename = fullfile(data_folder, data_files(k).name);
    
    % Read the data from the CSV file
    data = readmatrix(filename);
    
    % Create a time vector based on the sampling rate
    time = (0:length(data)-1)' / sampling_rate;
    
    % Apply the high-pass filter to the data
    filtered_data = filter(b, a, data);
    
    % Create a table with the current data, filtered data, and time
    current_table = table(time, data, filtered_data, 'VariableNames', {'Time', ['Data_' num2str(k)], ['Data_' num2str(k) '_Filtered']});
    
    % Append the current table to the combined data table
    if isempty(combined_data)
        combined_data = current_table;
    else
        combined_data = outerjoin(combined_data, current_table, 'Keys', 'Time', 'MergeKeys', true);
    end
end

% Plot all the original and filtered data in a single figure
for k = 2:2:width(combined_data)-1
    figure;
    hold on;
    plot(combined_data.Time, combined_data{:, k}, 'DisplayName', combined_data.Properties.VariableNames{k});
    plot(combined_data.Time, combined_data{:, k+1}, 'DisplayName', combined_data.Properties.VariableNames{k+1});
    hold off;
    xlabel('Time (seconds)');
    ylabel('Amplitude (mV)');
    title('Combined ECG Signals with Filtered Data');
    legend('show', 'Interpreter', 'none');
    grid on;
end

[pks, loc] = findpeaks(combined_data.Data_1_Filtered, 'MinPeakHeight', .2);
pksInverted = findpeaks(combined_data.Data_2_Filtered);

figure;
hold on;
plot(0:length(pks)-1, pks)
text(loc+.02,pks,num2str((1:numel(pks))'))
hold off;

figure;
plot(0:length(pksInverted)-1, pksInverted)


% Display the combined table
%disp(combined_data);
