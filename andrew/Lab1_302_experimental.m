%dataANDREWwristankleankle.acq
%Andrew Cunningham

%Define parameters
%fs = input('Input sampling frequency, in hz: '); 
duration = 60; 
%numSamples = fs * duration;

%get file
%filename = 'dataANDREWwristankleankle.txt';
%[filename, pathname] = uigetfile('*.txt', 'Select the ECG Data File');
%if isequal(filename, 0)
%    error('No file selected. Please select a valid file.');
%end
%filepath = fullfile(pathname, filename);
%data = readtable(filename, 'FileType', 'text');

% Get folder
foldername = uigetdir('', 'Select the folder containing ECG Data Files');
if isequal(foldername, 0)
    error('No folder selected. Please select a valid folder.');
end

% Get list of all .txt files in the folder
txtFiles = dir(fullfile(foldername, '*.txt'));

if isempty(txtFiles)
    error('No .txt files found in the selected folder.');
end

% Prompt user once for filter choice
applyFilter = input('Do you want to apply the filter to the ECG data (High-pass, 6Hz)? (y/n): ', 's');

% Initialize arrays to store results
average_BPMs = [];
std_BPMs = [];
sampling_frequencies = [];

for k = 1:length(txtFiles)
    filename = txtFiles(k).name;
    filepath = fullfile(foldername, filename);
    disp(['Processing file: ', filename]);

    freqMatch = regexp(filename, '\d+', 'match');
    if isempty(freqMatch)
        error(['No sampling frequency found in the file name: ', filename]);
    end
    fs = str2double(freqMatch{1});
    if isnan(fs)
        error(['Invalid sampling frequency extracted from the file name: ', filename]);
    end
    disp(['Extracted sampling frequency: ', num2str(fs), ' Hz']);  % Debugging print
    numSamples = fs * duration

    % Read the file using fopen and textscan
    %fid = fopen(filepath, 'r');
    %if fid == -1
    %    error(['Cannot open file: ', filename]);
    %end
    
    %data = textscan(fid, '%f');
    %fclose(fid);
    %ecgData = data{1};
    data = readtable(filepath, 'FileType', 'text');
    ecgData = table2array(data(:, 1));
    ecgData = ecgData(~isnan(ecgData));
    ecgData = table2array(data(:, 1));
ecgData = ecgData(~isnan(ecgData));

if length(ecgData) < numSamples
    error('Not enough data points in the file for the specified duration.');
end

% Create time vector
t = (0:numSamples-1) / fs;

% Extract the first 60 seconds of ECG data
ecgData = ecgData(1:numSamples);

% Plot ECG data
%figure;
%plot(t, ecgData);
%xlabel('Time (s)');
%ylabel('ECG Signal');
%title('ECG Data');
%grid on;

%figure;
%plot(t(1:1.5*fs), ecgData(1:1.5*fs)); % Plotting the first 1.5 seconds (1500 samples)
%xlabel('Time (s)');
%ylabel('ECG Signal');
%title('Single Waveform (First 1.5 Seconds)');
%grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END OF LAB 0%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply the filter to the ECG data
%applyFilter = input('Do you want tno apply the filter to the ECG data (High-pass, 6Hz)? (y/n): ', 's');

if strcmpi(applyFilter, 'y')
    %low-pass filter with a cutoff frequency of 6 Hz
    cutoff_freq = 6; 
    Wn = cutoff_freq / (fs / 2); 
    [b, a] = butter(4, Wn, 'high'); 
    % Apply the filter to the ECG data
    filtered_ecgData = filtfilt(b, a, ecgData);
    disp('Filter applied to the ECG data.');
else
    % Do not apply the filter, use original data
    filtered_ecgData = ecgData;
    disp('Filter not applied. Using original ECG data. Note that all peaks may not be identified.');
end

% Plot filtered ECG data
%figure;
%plot(t, filtered_ecgData);
%xlabel('Time (s)');
%ylabel('ECG Signal');
%title('Filtered ECG Data (Cutoff at 6 Hz)');
%grid on;

% Create a second figure for the filtered single waveform
%figure;
%plot(t(1:1.5*fs), filtered_ecgData(1:1.5*fs)); 
%xlabel('Time (s)');
%ylabel('ECG Signal');
%title('Filtered Single Waveform (First 1.5 Seconds, Cutoff at 6 Hz)');
%grid on;

% Find peaks in the filtered ECG data
[pks, locs] = findpeaks(filtered_ecgData, 'MinPeakHeight', 0.4, 'MinPeakProminence', 0.6);

% Calculate the vertical range and baseline
%verticalRange = max(filtered_ecgData) - min(filtered_ecgData);
%$baseline = min(filtered_ecgData) + 0.2 * verticalRange;
% Perform Fourier transform to estimate baseline
fft_ecg = fft(filtered_ecgData);
freqs = (0:numSamples-1) * (fs / numSamples);

% Set high frequency components to zero to isolate low-frequency baseline
cutoff_freq_baseline = 0.5; % 0.5 Hz as an example cutoff for baseline
fft_ecg(abs(freqs) > cutoff_freq_baseline) = 0;

% Perform inverse Fourier transform to get the baseline
baseline_ecg = ifft(fft_ecg, 'symmetric');

% Calculate average amplitude of the peaks relative to the baseline
amplitudes = pks - baseline_ecg(locs);
average_amplitude = mean(amplitudes);

% Plot the filtered ECG data with peaks highlighted
%figure;
%plot(t, filtered_ecgData);
%hold on;
%plot(t(locs), pks, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
%xlabel('Time (s)');
%ylabel('ECG Signal');
%title('ECG Data with Peaks');
%legend('ECG', 'Peaks');
%grid on;

% Create a second figure for the filtered single waveform with peaks highlighted
%figure;
%plot(t(1:1.5*fs), filtered_ecgData(1:1.5*fs));
%hold on;
%single_waveform_locs = locs(locs <= 1.5*fs);
%single_waveform_pks = pks(locs <= 1.5*fs);
%plot(t(single_waveform_locs), single_waveform_pks, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
%xlabel('Time (s)');
%ylabel('ECG Signal');
%title('Single Waveform with Peaks (First 1.5 Seconds)');
%legend('ECG', 'Peaks');
%grid on;

%figure;
%plot(t(1:20*fs), filtered_ecgData(1:20*fs));
%hold on;
%single_waveform_locs = locs(locs <= 20*fs);
%single_waveform_pks = pks(locs <= 20*fs);
%plot(t(single_waveform_locs), single_waveform_pks, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
%xlabel('Time (s)');
%ylabel('ECG Signal');
%title('Single Waveform with Peaks (First 20 Seconds)');
%legend('ECG', 'Peaks');
%grid on;

% Calculate R-R intervals
RR_intervals = diff(t(locs));
% Calculate BPM for each R-R interval
BPM = 60 ./ RR_intervals;
% Calculate average BPM and standard deviation 
average_BPM = mean(BPM);
std_BPM = std(BPM);
% Number of heartbeats sampled
num_heartbeats = length(RR_intervals);

% Display R-R intervals
disp('R-R Intervals (seconds):');
disp(RR_intervals);
disp(['Average BPM: ', num2str(average_BPM)]);
disp(['Standard Deviation of BPM: (+/-)', num2str(std_BPM)]);
disp(['Number of Heartbeats Sampled: ', num2str(num_heartbeats)]);
disp(['Average Amplitude of Peaks: ', num2str(average_amplitude)]);

% Store results
    average_BPMs(end+1) = average_BPM;
    std_BPMs(end+1) = std_BPM;
    sampling_frequencies(end+1) = fs;

%disp(['shut up no one cares about your warning matlab'])

%Plot the filtered ECG data with peaks and baseline highlighted
figure;
plot(t, filtered_ecgData);
hold on;
plot(t(locs), pks, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
plot(t, baseline_ecg, '--g', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('ECG Signal');
title(['ECG Data with Peaks and Baseline: ', ' (Sampling Frequency: ', num2str(fs), ' Hz)']);
legend('ECG', 'Peaks', 'Baseline');
grid on;
% Add text box with statistics
annotation('textbox', [0.15, 0.7, 0.2, 0.2], 'String', {['Average BPM: ', num2str(average_BPM)], ...
    ['Standard Deviation of BPM: ', num2str(std_BPM)], ...
    ['Number of Heartbeats Sampled: ', num2str(num_heartbeats)], ...
    ['Average Amplitude of Peaks (millivolts): ', num2str(average_amplitude)]}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white');
hold off

end

% Plot standard deviation versus average BPM
figure;
scatter(average_BPMs, std_BPMs, 'filled');
xlabel('Average BPM');
ylabel('Standard Deviation of BPM');
title('Standard Deviation vs. Average BPM');
grid on;

% Annotate points with file names
for i = 1:length(sampling_frequencies)
    text(average_BPMs(i), std_BPMs(i), num2str(sampling_frequencies(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end





