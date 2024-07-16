% Get the full path of the folder containing this script
script_folder = fileparts(mfilename('fullpath'));

% Define the path to the data folder relative to the script folder
data_folder = fullfile(script_folder, '../lab_0_data/converted_csvs/');

% Get a list of all text files in the data folder
data_files = dir(fullfile(data_folder, '*.txt'));

% Loop through each file and process it
for k = 1:length(data_files)
    % Get the full path of the current file
    filename = fullfile(data_folder, data_files(k).name);
    
    % Open the file for reading
    fileID = fopen(filename, 'r');
    
    % Read the first line (discard)
    fgetl(fileID);
    
    % Read the second line to extract the sampling rate
    sampling_rate_line = fgetl(fileID);
    % Extract the numeric part of the sampling rate (assume it's in msec/sample)
    sampling_rate_msec = sscanf(sampling_rate_line, '%f');
    % Convert sampling rate to Hz (samples per second)
    sampling_rate = 1000 / sampling_rate_msec;
    
    % Read and discard the next 5 lines of the header
    for i = 1:5
        fgetl(fileID);
    end
    
    % Read the data from the file
    data = fscanf(fileID, '%f');
    
    % Close the file
    fclose(fileID);
    
    % Create a time vector based on the sampling rate
    time = (0:length(data)-1) / sampling_rate;
    
    % Plot the data
    figure;
    plot(time, data);
    xlabel('Time (seconds)');
    ylabel('Amplitude (mV)');
    title(['ECG Signal - ' data_files(k).name]);
    grid on;
end
