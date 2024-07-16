% Function to find peaks and calculate BPM
function [peaks, bpm] = findPeaks(data, time, frequency)
    % Find peaks in the data
    [peaks, locs] = findpeaks(data, 'MinPeakHeight', mean(data) + std(data));
    
    % Calculate the time between peaks
    peakIntervals = diff(time(locs));
    
    % Calculate BPM
    avgPeakInterval = mean(peakIntervals);
    bpm = 60 / avgPeakInterval;
end