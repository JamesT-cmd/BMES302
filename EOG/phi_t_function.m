function phi_t = phi_t_function(V, V1, V2)
    % PHI_T_FUNCTION Calculate the phi(t) values for the given EOG data
    %
    % Inputs:
    %   V - EOG signal (vector)
    %   V1 - Average of closest positive peaks
    %   V2 - Average of closest negative peaks
    %
    % Output:
    %   phi_t - Calculated phi(t) values (vector)
    %
    % Calculated  V1 and V2 for initial lab tests
    % V1: .6070         V2: -.5668

    if nargin < 2
        V1 = 0.6070;
    end
    if nargin < 3
        V2 = -.5668;
    end


    phi_t = 90 / (V1 - V2) * (V - (V1 + V2) / 2);
end
