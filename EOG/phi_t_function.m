function phi_t = phi_t_function(V, V1_or_calibration, V2)
    % PHI_T_FUNCTION Calculate the phi(t) values for the given EOG data
    %
    % Inputs:
    %   V - EOG signal (vector)
    %   V1_or_calibration - Average of closest positive peaks or calibration set number (1 or 2)
    %   V2 - Average of closest negative peaks (optional, used only if V1 is supplied)
    %
    % Output:
    %   phi_t - Calculated phi(t) values (vector)
    %
    % Calculated V1 and V2 for initial lab tests
    % V1: .6070         V2: -.5668
    %
    % Calculated V1 and V2 for second calibrations
    % V1: .5814         V2: -.5233

    if nargin == 2  % If only V1_or_calibration is supplied, use it as calibration
        calibration = V1_or_calibration;
        if calibration == 1
            V1 = 0.6070;
            V2 = -0.5668;
        elseif calibration == 2
            V1 = 0.5814;
            V2 = -0.5233;
        else
            error('Invalid calibration number. Use 1 or 2.');
        end
    elseif nargin == 3  % If V1 and V2 are supplied
        V1 = V1_or_calibration;
        V2 = V2;
    else
        error('Invalid number of arguments. Provide either (V, calibration) or (V, V1, V2).');
    end

    phi_t = 90 / (V1 - V2) * (V - (V1 + V2) / 2);
end
