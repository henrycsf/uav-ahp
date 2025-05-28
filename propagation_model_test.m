%% Simulation Parameters
N_users = 400;                  % Number of users
area_size = 2000;               % Dimensions of the area (2000x2000 m)
h = 100;                        % UAV altitude in meters

% UAV positions (placed at the centers of the four quadrants)
UAV_positions = [ area_size/4, area_size/4;       % UAV 1: bottom left
                  3*area_size/4, area_size/4;     % UAV 2: bottom right
                  area_size/4, 3*area_size/4;     % UAV 3: top left
                  3*area_size/4, 3*area_size/4];  % UAV 4: top right
N_UAV = size(UAV_positions,1);

% Propagation Parameters
f = 3.5e9;                      % Frequency in Hz
velc = 299792458;               % Speed of light (m/s)
ZetaLOS = 1;                    % LOS loss factor (dB)
ZetaNLOS = 20;                  % NLOS loss factor (dB)
alpha = 4.88;                   % Model parameter alpha
beta = 0.83;                    % Model parameter beta

% Transmit parameters
Pt = 35;                        % Transmit power in dBm
Gt = 3;                         % Transmitter gain in dBi
Gr = 0;                         % Receiver gain in dBi

% Noise calculation (using a bandwidth of 50 MHz)
B = 50e6;                       % Bandwidth in Hz
q = -174 + 10*log10(B);          % Noise power in dBm for bandwidth B
noise_lin = 10^((q-30)/10);      % Convert noise power to linear scale (W)

%% Generate User Locations
% Uniformly distributed in the area [0,2000] x [0,2000]
users = area_size * rand(N_users, 2);

%% Preallocate arrays for received powers (rows: users, columns: UAVs)
Pr_dBm = zeros(N_users, N_UAV);
Pr_lin = zeros(N_users, N_UAV);

%% Calculate Received Power from Each UAV
for b = 1:N_UAV
    % Get UAV b position
    UAV_pos = UAV_positions(b,:);
    
    % Horizontal distance R from UAV to each user
    R = sqrt((users(:,1) - UAV_pos(1)).^2 + (users(:,2) - UAV_pos(2)).^2);
    % 3D distance including altitude
    D = sqrt(R.^2 + h^2);
    
    % Elevation angle theta (in radians)
    % Using atan2 to handle R==0 properly
    theta = atan2(h, R);
    theta(R==0) = pi/2;
    
    % Additional loss factor Z (angle converted to degrees)
    Z = alpha * exp(-beta * ((theta*180/pi) - alpha));
    
    % Path Loss (PL) in dB: free-space term + additional loss term
    PL = 20*log10((4*pi*f.*D)/velc) + ((ZetaLOS + Z.*ZetaNLOS) ./ (1 + Z));
    
    % Received power in dBm for UAV b
    Pr_dBm(:,b) = Pt - PL + Gt + Gr;
    
    % Convert received power to linear scale (W)
    Pr_lin(:,b) = 10.^((Pr_dBm(:,b)-30)/10);
end

%% SINR Calculation Considering Interference from the 3 Other UAVs
SINR_lin = zeros(N_users, N_UAV);
for b = 1:N_UAV
    % Sum the interference from all other UAVs
    interference = sum(Pr_lin(:,[1:b-1, b+1:end]), 2);
    % Calculate SINR for UAV b: signal power divided by (noise + interference)
    SINR_lin(:,b) = Pr_lin(:,b) ./ (noise_lin + interference);
end

% For each user, choose the UAV that provides the maximum SINR
[max_SINR_lin, best_UAV] = max(SINR_lin, [], 2);
max_SINR_dB = 10*log10(max_SINR_lin);

%% Plotting the SINR Colormap
figure;
scatter(users(:,1), users(:,2), 50, max_SINR_dB, 'filled');
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
title('User SINR (dB) with 4 Interfering UAVs');
axis([0 area_size 0 area_size]);
axis equal;
grid on;
