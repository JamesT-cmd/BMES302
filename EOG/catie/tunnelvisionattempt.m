% MATLAB script to read and plot baseline and tunnel vision EOG data from text files

% Define the sampling interval (in milliseconds)
samplingInterval = 1; % Adjust this value based on your actual sampling rate

% Read baseline EOG data
baselineEOG = readmatrix('baseline2.txt');
timeBaseline = (0:length(baselineEOG)-1) * samplingInterval;

% Truncate baseline data at 4000 ms
endTime = 4000; % End time in milliseconds
endIndexBaseline = find(timeBaseline <= endTime, 1, 'last');
timeBaseline = timeBaseline(1:endIndexBaseline);
baselineEOG = baselineEOG(1:endIndexBaseline);


% Read tunnel vision EOG data
tunnelVisionEOG = readmatrix('tunnelvision.txt');
timeTunnelVision = (0:length(tunnelVisionEOG)-1) * samplingInterval;

% Truncate tunnel vision data at 4000 ms
endIndexTunnelVision = find(timeTunnelVision <= endTime, 1, 'last');
timeTunnelVision = timeTunnelVision(1:endIndexTunnelVision);
tunnelVisionEOG = tunnelVisionEOG(1:endIndexTunnelVision);

numSubplots = 2;

%plot
figure;
subplot(numSubplots, 1, 1);
plot(timeBaseline, baselineEOG, 'b-', 'LineWidth', 1.5);
title('Baseline EOG');
xlabel('Time (ms)');
ylabel('EOG Signal (uV)');
grid on;

subplot(numSubplots, 1, 2);
plot(timeTunnelVision, tunnelVisionEOG, 'r-', 'LineWidth', 1.5);
title('Tunnel Vision EOG');
xlabel('Time (ms)');
ylabel('EOG Signal (uV)');
grid on;
% MATLAB script to plot baseline EOG data against tunnel vision EOG data

% Construct the relative path to the data file
%relative_path = ['..', filesep, 'lab_3_data', filesep, 'baseline2.txt'];

% Load the data from the relative path
%data = load(relative_path);

%baselineEOG = 'baseline2.txt';

%tunnelVisionEOG = 'tunnelvision.txt';

% Time vector (assuming sampling frequency is the same for both datasets)
%time = 1:length(baselineEOG);  % Replace with appropriate time vector if needed

% Plot the data
%figure;
%plot(time, baselineEOG, 'b-', 'LineWidth', 1.5); % Plot baseline EOG in blue
%hold on;
%plot(time, tunnelVisionEOG, 'r-', 'LineWidth', 1.5); % Plot tunnel vision EOG in red
%hold off;

% Add titles and labels
%title('Baseline EOG vs Tunnel Vision EOG');
%xlabel('Time (s)');
%ylabel('EOG Signal (uV)');
%legend('Baseline EOG', 'Tunnel Vision EOG');
%grid on;