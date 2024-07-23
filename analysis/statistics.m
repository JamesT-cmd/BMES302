% Extract data for ANDREW2 from resultsTable0
andrew2_0 = resultsTable0(strcmp(resultsTable0.Name, 'ANDREW2'), :);

% Display the extracted data for verification
disp('Data for ANDREW2 from resultsTable0:');
disp(andrew2_0);

% Initialize arrays to store p-values and h-values from t-tests
p_heart_rate = zeros(height(resultsTable2), 1);
h_heart_rate = zeros(height(resultsTable2), 1);
p_peak_size = zeros(height(resultsTable2), 1);
h_peak_size = zeros(height(resultsTable2), 1);

% Perform statistical tests for each row in resultsTable2 using ANDREW2 in resultsTable0
for i = 1:height(resultsTable2)
    % Extract data for the current row in resultsTable2
    heart_rate_2 = resultsTable2.BPM(i);
    std_heart_rate_2 = resultsTable2.StdPeakInterval(i);
    peak_size_2 = resultsTable2.AveragePeakHeight(i);
    std_peak_size_2 = resultsTable2.StdPeakHeight(i);
    total_peaks_2 = resultsTable2.TotalPeaks(i);
    
    % Extract data for ANDREW2 from resultsTable0
    heart_rate_0 = andrew2_0.AverageHeartRate;
    std_heart_rate_0 = andrew2_0.StdHeartRate;
    peak_size_0 = andrew2_0.AveragePeakMagnitude;
    std_peak_size_0 = andrew2_0.StdPeakHeight;
    total_peaks_0 = andrew2_0.TotalPeaks;

    % Calculate the pooled standard deviation for heart rate
    pooled_std_heart_rate = sqrt(((total_peaks_0 - 1) * std_heart_rate_0^2 + (total_peaks_2 - 1) * std_heart_rate_2^2) / (total_peaks_0 + total_peaks_2 - 2));
    % Calculate the t-statistic for heart rate
    t_heart_rate = abs(heart_rate_0 - heart_rate_2) / (pooled_std_heart_rate * sqrt(1/total_peaks_0 + 1/total_peaks_2));
    % Calculate the degrees of freedom for heart rate
    df_heart_rate = total_peaks_0 + total_peaks_2 - 2;
    % Calculate the p-value for heart rate
    p_heart_rate(i) = 2 * tcdf(-abs(t_heart_rate), df_heart_rate);
    % Determine if the result is significant
    h_heart_rate(i) = p_heart_rate(i) < 0.05;
    
    % Calculate the pooled standard deviation for peak size
    pooled_std_peak_size = sqrt(((total_peaks_0 - 1) * std_peak_size_0^2 + (total_peaks_2 - 1) * std_peak_size_2^2) / (total_peaks_0 + total_peaks_2 - 2));
    % Calculate the t-statistic for peak size
    t_peak_size = abs(peak_size_0 - peak_size_2) / (pooled_std_peak_size * sqrt(1/total_peaks_0 + 1/total_peaks_2));
    % Calculate the degrees of freedom for peak size
    df_peak_size = total_peaks_0 + total_peaks_2 - 2;
    % Calculate the p-value for peak size
    p_peak_size(i) = 2 * tcdf(-abs(t_peak_size), df_peak_size);
    % Determine if the result is significant
    h_peak_size(i) = p_peak_size(i) < 0.05;
end

% Create a table to store the t-test results
ttestResults = table(resultsTable2.Frequency, resultsTable2.BPM, resultsTable2.AveragePeakHeight, ...
    h_heart_rate, p_heart_rate, h_peak_size, p_peak_size, ...
    'VariableNames', {'Frequency', 'BPM', 'AveragePeakHeight', 'HeartRate_h', 'HeartRate_p', 'PeakSize_h', 'PeakSize_p'});

% Display the t-test results
disp('t-test Results for Each Row in resultsTable2 Compared to ANDREW2 from resultsTable0:');
disp(ttestResults);

% Save the t-test results to the workspace
assignin('base', 'ttestResults', ttestResults);

% Extract data for the 1000 Hz frequency from resultsTable2
data_1000Hz = resultsTable2(resultsTable2.Frequency == 1000, :);

% Initialize arrays to store p-values and h-values from t-tests
p_heart_rate_1000Hz = zeros(height(resultsTable2), 1);
h_heart_rate_1000Hz = zeros(height(resultsTable2), 1);
p_peak_size_1000Hz = zeros(height(resultsTable2), 1);
h_peak_size_1000Hz = zeros(height(resultsTable2), 1);

% Perform statistical tests for each row in resultsTable2 compared to 1000 Hz frequency data
for i = 1:height(resultsTable2)
    if resultsTable2.Frequency(i) == 1000
        continue; % Skip comparison with itself
    end
    
    % Extract data for the current row in resultsTable2
    heart_rate_i = resultsTable2.BPM(i);
    std_heart_rate_i = resultsTable2.StdPeakInterval(i);
    peak_size_i = resultsTable2.AveragePeakHeight(i);
    std_peak_size_i = resultsTable2.StdPeakHeight(i);
    total_peaks_i = resultsTable2.TotalPeaks(i);

    % Perform t-tests using the number of peaks as the sample size
    % Calculate the pooled standard deviation for heart rate
    pooled_std_heart_rate = sqrt(((data_1000Hz.TotalPeaks - 1) * data_1000Hz.StdPeakInterval^2 + (total_peaks_i - 1) * std_heart_rate_i^2) / (data_1000Hz.TotalPeaks + total_peaks_i - 2));
    % Calculate the t-statistic for heart rate
    t_heart_rate = abs(data_1000Hz.BPM - heart_rate_i) / (pooled_std_heart_rate * sqrt(1/data_1000Hz.TotalPeaks + 1/total_peaks_i));
    % Calculate the degrees of freedom for heart rate
    df_heart_rate = data_1000Hz.TotalPeaks + total_peaks_i - 2;
    % Calculate the p-value for heart rate
    p_heart_rate_1000Hz(i) = 2 * tcdf(-abs(t_heart_rate), df_heart_rate);
    % Determine if the result is significant
    h_heart_rate_1000Hz(i) = p_heart_rate_1000Hz(i) < 0.05;
    
    % Calculate the pooled standard deviation for peak size
    pooled_std_peak_size = sqrt(((data_1000Hz.TotalPeaks - 1) * data_1000Hz.StdPeakHeight^2 + (total_peaks_i - 1) * std_peak_size_i^2) / (data_1000Hz.TotalPeaks + total_peaks_i - 2));
    % Calculate the t-statistic for peak size
    t_peak_size = abs(data_1000Hz.AveragePeakHeight - peak_size_i) / (pooled_std_peak_size * sqrt(1/data_1000Hz.TotalPeaks + 1/total_peaks_i));
    % Calculate the degrees of freedom for peak size
    df_peak_size = data_1000Hz.TotalPeaks + total_peaks_i - 2;
    % Calculate the p-value for peak size
    p_peak_size_1000Hz(i) = 2 * tcdf(-abs(t_peak_size), df_peak_size);
    % Determine if the result is significant
    h_peak_size_1000Hz(i) = p_peak_size_1000Hz(i) < 0.05;
end

% Create a table to store the t-test results for 1000 Hz comparisons
ttestResults_1000Hz = table(resultsTable2.Frequency, resultsTable2.BPM, resultsTable2.AveragePeakHeight, ...
    h_heart_rate_1000Hz, p_heart_rate_1000Hz, h_peak_size_1000Hz, p_peak_size_1000Hz, ...
    'VariableNames', {'Frequency', 'BPM', 'AveragePeakHeight', 'HeartRate_h_1000Hz', 'HeartRate_p_1000Hz', 'PeakSize_h_1000Hz', 'PeakSize_p_1000Hz'});

% Display the t-test results for 1000 Hz comparisons
disp('t-test Results for 1000 Hz Frequency Compared to Other Rows in resultsTable2:');
disp(ttestResults_1000Hz);

% Save the t-test results to the workspace
assignin('base', 'ttestResults_1000Hz', ttestResults_1000Hz);
