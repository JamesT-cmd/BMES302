% Define the relative path to the data folder
data_folder = ['..', filesep, 'lab_3_data'];

% Get a list of all .txt files in the folder
files = dir(fullfile(data_folder, '*.txt'));

% Sampling rate
sampling_rate = 200;

% Loop through each file
for k = 1:length(files)
    % Construct the full file path
    file_path = fullfile(files(k).folder, files(k).name);
    
    % Load the data from the file
    data = load(file_path);
    
    % Calculate the time vector
    time = (0:length(data)-1) / sampling_rate;
    
    % Compute phi(t) using the function
    phi_t = phi_t_function(data);
    
    % Plot phi(t) over time
    figure;
    plot(time, phi_t, 'm', 'DisplayName', '\phi(t)');
    hold on;
    yline(45, 'r--', 'DisplayName', '45 Degrees');
    yline(-45, 'r--', 'DisplayName', '-45 Degrees');
    xlabel('Time (s)');
    ylabel('\phi(t)');
    legend show;
    title(['\phi(t) Calculation for ', files(k).name]);
    grid on;
end
