% Define the relative path to the data folders
data_folder = fullfile('..', 'lab_3_data');

% Define the directory paths for the calibration sets
calibration1_dir = fullfile(data_folder, 'calibration1');
calibration2_dir = fullfile(data_folder, 'calibration2');

% Get a list of all .txt files in each calibration folder
calibration1_files = dir(fullfile(calibration1_dir, '*.txt'));
calibration2_files = dir(fullfile(calibration2_dir, '*.txt'));

% Sampling rate
sampling_rate = 200;

% Loop through each file in the calibration1 directory
for k = 1:length(calibration1_files)
    % Construct the full file path
    file_path = fullfile(calibration1_files(k).folder, calibration1_files(k).name);
    
    % Load the data from the file
    data = load(file_path);
    
    % Process and plot the data
    process_data(data, sampling_rate, 1, calibration1_files(k).name);
end

% Loop through each file in the calibration2 directory
for k = 1:length(calibration2_files)
    % Construct the full file path
    file_path = fullfile(calibration2_files(k).folder, calibration2_files(k).name);
    
    % Load the data from the file
    data = load(file_path);
    
    % Process and plot the data
    process_data(data, sampling_rate, 2, calibration2_files(k).name);
end
