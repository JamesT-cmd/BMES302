function process_data(data, sampling_rate, calibration, file_name)
    % Calculate the time vector
    time = (0:length(data)-1) / sampling_rate;
    
    % Compute phi(t) using the function
    phi_t = phi_t_function(data, calibration);
    
    % Plot phi(t) over time
    figure;
    plot(time, phi_t, 'm', 'DisplayName', '\phi(t)');
    hold on;
    yline(45, 'r--', 'DisplayName', '45 Degrees');
    yline(-45, 'r--', 'DisplayName', '-45 Degrees');
    xlabel('Time (s)');
    ylabel('\phi(t)');
    legend show;
    title(['\phi(t) Calculation for ', file_name]);
    grid on;
end