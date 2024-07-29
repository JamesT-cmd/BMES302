% Construct the relative path to the data file
relative_path = ['..', filesep, 'lab_3_data', filesep, 'calibrationSetLeftFirst.txt'];

% Load the data from the relative path
data = load(relative_path);

mu = mean(data);
sigma = std(data);

% Define thresholds for peak detection
threshold = mu + 2 * sigma;

% Find local maxima (positive spikes) above the positive_threshold
[positive_peaks, pos_locs] = findpeaks(data, 'MinPeakHeight', threshold);

% Find local minima (negative spikes) below the negative_threshold
[negative_peaks, neg_locs] = findpeaks(-data, 'MinPeakHeight', threshold);

% Convert sample indices to time (in seconds)
time = (0:length(data)-1) / 200;
pos_time = pos_locs / 200;
neg_time = neg_locs / 200;

% Plot the raw data with all detected peaks
figure;
plot(time, data, 'b', 'DisplayName', 'Raw Data');
hold on;
plot(pos_time, data(pos_locs), 'g^', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Positive Peaks');
plot(neg_time, data(neg_locs), 'rv', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Negative Peaks');
xlabel('Time (s)');
ylabel('EOG Signal');
legend show;
title('Detected Peaks in EOG Data');
grid on;

% Get the leftmost side of the positive spikes
left_side_positive_spikes = [];
for i = 1:length(pos_locs)
    peak_index = pos_locs(i);
    % Move left until the value falls below the threshold
    while peak_index > 1 && data(peak_index - 1) > threshold
        peak_index = peak_index - 1;
    end
    left_side_positive_spikes = [left_side_positive_spikes; peak_index];
end

% Get the rightmost side of the negative spikes
right_side_negative_spikes = [];
for i = 1:length(neg_locs)
    peak_index = neg_locs(i);
    % Move right until the value falls below the threshold
    while peak_index < length(data) && data(peak_index + 1) < -threshold
        peak_index = peak_index + 1;
    end
    right_side_negative_spikes = [right_side_negative_spikes; peak_index];
end

% Convert spike indices to time
left_side_positive_time = left_side_positive_spikes / 200;
right_side_negative_time = right_side_negative_spikes / 200;

% Find the closest positive and negative peaks
closest_positive_peaks = [];
closest_negative_peaks = [];

for i = 1:length(left_side_positive_spikes)
    [~, idx] = min(abs(pos_locs - left_side_positive_spikes(i)));
    closest_positive_peaks = [closest_positive_peaks; pos_locs(idx)];
end

for i = 1:length(right_side_negative_spikes)
    [~, idx] = min(abs(neg_locs - right_side_negative_spikes(i)));
    closest_negative_peaks = [closest_negative_peaks; neg_locs(idx)];
end

% Convert closest peak indices to time
closest_positive_time = closest_positive_peaks / 200;
closest_negative_time = closest_negative_peaks / 200;

% Calculate V1 and V2
V1 = mean(data(closest_positive_peaks));
V2 = mean(data(closest_negative_peaks));

% Compute phi(t) using the function
phi_t = phi_t_function(data);

% Plot the results with the selected positive and negative points
figure;
plot(time, data, 'b', 'DisplayName', 'Raw Data');
hold on;
plot(closest_positive_time, data(closest_positive_peaks), 'go', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Closest Positive Peaks');
plot(closest_negative_time, data(closest_negative_peaks), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Closest Negative Peaks');
xlabel('Time (s)');
ylabel('EOG Signal');
legend show;
title('High End Drift Detection in EOG Data');
grid on;

% Plot phi(t)
figure;
plot(time, phi_t, 'm', 'DisplayName', '\phi(t)');
hold on;
yline(45, 'r--', 'DisplayName', '45 Degrees', 'LineWidth', 2);
yline(-45, 'r--', 'DisplayName', '-45 Degrees', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('\phi(t)');
legend show;
title('\phi(t) Calculation');
grid on;

% Display V1, V2, and phi(t) values
disp('V1 (Average of Closest Positive Peaks):');
disp(V1);
disp('V2 (Average of Closest Negative Peaks):');
disp(V2);
disp('phi(t):');
disp(phi_t);
