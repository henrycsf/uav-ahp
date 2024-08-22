function [Zcov,Zcap,R_cov,R_cap,Pt] = numerology(h,X,Y,N_users,delta_cov,DataRate, alpha, beta)

warning('off','all');

W = N_users;
csi = 0.5;

A = 1.95;
Area = X * Y * 10^-6;
User_Density = N_users / (Area);

%% 5G Numerology

Mu = [0 1 2];
Latency = [1*10^-3 0.5*10^-3 0.25*10^-3];
delta_f = [15*10^3 30*10^3 60*10^3];
B = [50*10^6 20*10^6 10*10^6];
vlayers = [4 4 2];
Qm = [8 8 6];

%% Maximum Transmission Rate

factor = 1;
Rmax = 948/1024;
OH = 0.14;
Ts = 10^-3 ./ (14*2.^Mu);

N_RB = [270 51 11];

V = vlayers .* Qm .* factor .* Rmax;

TxRate = 10^-6 * (V .* N_RB .* 12 .* (1-OH))./Ts;

p_connected = 0.2;
p_nonactive = 0.15;

p_active = p_connected + p_nonactive;

Wx_thr = cap_modelling(W, csi, TxRate, p_active, DataRate);

N_SC = 12;
Nx_RB = N_RB;
Nthr_RB = sum(N_RB);
Nx = N_RB .* N_SC;
Nx_CET = floor(0.95 .* Nx_RB);

delta_r = sqrt(X^2 + Y^2)/20;
delta_l = 0.05;

site_cap = 3 * sum(Wx_thr);

syms R;

f = 3.5;
velc = 299792458;
ZetaLOS = 1;
ZetaNLOS = 20;
Gt = 3;
Gr = 0;

% h = 50;
% h_UE = 1.5;
% 
% d_BP = 4 * (h_UE - 1) * (h-1) * f * 10^9/ velc;

pen_loss = 3;
S_UE = -90;
sigma = 4;

IL = (-10) * log10(1 - delta_cov);

L = IL + pen_loss + S_UE + sigma;

Pmax = 40;

Rcov_x = zeros(1,3);

Ex = zeros(1,length(Nx_CET));
Ex_CET = zeros(1,length(Nx_CET));
PL_max = zeros(1,length(Ex_CET));

for Ptx = 10:1:Pmax
    for x = 1:length(TxRate)
            
            Ex(x) = Ptx - 10.*log10(Nx(x));

            Ex_CET(x) = Ex(x) + 10 * log10(Nx_CET(x) * N_SC);

            PL_max(x) = Ex_CET(x) - L;

            theta = atan(h./R);

            Z = (alpha*exp(-beta*((180/pi).*theta - alpha))); 
            
            PL = 20*log10((4*pi*(f*10^9)*R./velc))+((ZetaLOS+Z.*ZetaNLOS)./(1+Z));
            
            Rcov_solve = solve(PL == PL_max(x), R) * 10^-3;
    
            Rcov_x(x) = cast(Rcov_solve(1), 'double');
            
            Rcov_x(x) = abs(Rcov_x(x));
    end

    Rcov = min(Rcov_x);

    for delta_cap = 0:0.01:1

        Ac = (site_cap * delta_cap / User_Density);
        Rcap = (sqrt(Ac/A));

        delta_R = 10^3 * abs(Rcov - Rcap);
        delta_L = abs(delta_cov - delta_cap);

        if (delta_R<=delta_r && delta_L <= delta_l)
            Zcov = ceil(Area / (A*Rcov^2));
            Zcap = ceil(Area / Ac);
            R_cov = Rcov;
            R_cap = Rcap;
            Pt = Ptx;
        end
    end
end

end

function Wx_thr = cap_modelling(W, csi, TxRate, p_active, DataRate)

    Wx_thr = zeros(1,length(TxRate));

    Wx_max = floor((TxRate .* 10^6) ./ DataRate);

    for i = 1:length(Wx_max)

        for j = Wx_max(i):1:W

            exp = binocdf(Wx_max(i),j,p_active);

            if exp >= csi-0.1 && exp <= csi
            
                Wx_thr(i) = j;
                break;
            end
        end

    end

end