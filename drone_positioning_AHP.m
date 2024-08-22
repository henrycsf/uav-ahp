function [UAV, Best, AHP_time, Assoc, Avg_SINR_final, Avg_TP_final] = drone_positioning_AHP(CaseSelect, Pt, N_users, X, Y, h, DataRate, ...
    Users, Demand, Zcov, Zcap, Rcov, Rcap, nmi, alpha, beta, plots)

tic

R_limit = min(Rcov, Rcap)*10^3;

UAV_limit = max(Zcov, Zcap);

%% Initial Solutions

ab=[X Y];

for ik = 1:1:UAV_limit
    
    Radius(ik,1) = R_limit;
    
end

Center = [];
cnst=true;
vis=false;

Step = 0.3*[X/max([X Y]) * R_limit, Y/(max([X Y])) * R_limit];

if R_limit < X && R_limit < Y
    dist = (2*R_limit/sqrt(X^2+Y^2));
else
    dist = 0;
end

while size(Center,1) ~= UAV_limit
[Center] = ([(1-(1/UAV_limit))*X (1-(1/UAV_limit))*Y] .* random_min_spacing(UAV_limit, dist, vis)) + Step;
end

N_UAV = length(Radius);

switch CaseSelect
    case 1
    C = [5; 7; 3; 5; 1/3; 1/5];
    case 2
    C = [1/5; 3; 1/3; 7; 5; 3];
    case 3
    C = [3; 3; 1/7; 3; 1/9; 1/9];

end

for i = 1:1:N_UAV

    Points_Linear(:,:,i) = [Center(i,1) Center(i,2);
        Center(i,1) Center(i,2)+Step(2);
        Center(i,1)-Step(1) Center(i,2)+Step(2);
        Center(i,1)-Step(1) Center(i,2);
        Center(i,1)-Step(1) Center(i,2)-Step(2);
        Center(i,1) Center(i,2)-Step(2);
        Center(i,1)+Step(1) Center(i,2)-Step(2);
        Center(i,1)+Step(1) Center(i,2);
        Center(i,1)+Step(1) Center(i,2)+Step(2)
    ];
    
    for j = 1:size(Points_Linear,1)
        
        UAV(j,:,i) = [Points_Linear(j,1,i) Points_Linear(j,2,i)];
        
    end

Pos_points = size(Points_Linear,1);

User = zeros(1,1,N_UAV);

I = ones(1,1,N_UAV);

Sol = [User I];

Best = [User I];

%% Calculating Network Demands

    f = 3.5 * 10^9;
    velc = 299792458;
    ZetaLOS = 1;
    ZetaNLOS = 20;
    BW = [50*10^6 20*10^6 10*10^6];
    q = -174 + 10*log10(BW);

    Gt = 3;
    Gr = 0;

    for j = 1:size(Users,1)

        for k = 1:1:Pos_points
            D(j,k,i) = sqrt(((Users(j,1)-UAV(k,1,i)).^2) + (Users(j,2)-UAV(k,2,i)).^2 + h.^2);
            R(j,k,i) = sqrt(abs(((Users(j,1)-UAV(k,1,i)).^2) + (Users(j,2)-UAV(k,2,i)).^2));
            Demand_Density(j,k,i) = Demand(j)/(D(j,k,i).^2);
        end

        switch Demand(j)
             case DataRate(1)
                 User_Service(j) = 1;
                 
             case DataRate(2)
                 User_Service(j) = 2;
                 
             case DataRate(3)
                 User_Service(j) = 3;

             otherwise
                 error("Out of service");
        end

    end

    theta = atan(h./R);

    Z = (alpha*exp(-beta*((180/pi).*theta - alpha)));

    PL = 20*log10((4*pi*f*D./velc))+((ZetaLOS+Z.*ZetaNLOS)./(1+Z));

    Pr_user = Pt - PL + Gt + Gr;

    Pr_Linear = (10.^((Pr_user-30)./10));

end

%% Iterative process begins

n = 0;

tol = 0;

while (n < nmi)

    n = n+1;

    Assoc_Matrix = zeros(1,Pos_points,N_UAV);

for i = 1:1:N_UAV

    for b = 1:1:Pos_points

            for j = 1:1:N_users

            inter = 0;

                for m = 1:1:N_UAV

                    if (i ~= m)
    
                        inter = inter + Pr_Linear(j,Best(1,2,i),m);
    
                    end
                end

            Interference_lin(j,b,i) = inter;

            Interference(j,b,i) = 10*log10(inter);

            switch Demand(j)

                case DataRate(1)
                if R(j,b,i) > R_limit
                    SINR_lin(j,b,i) = 0;
                    TP(j,b,i) = 0;
                else
                    SINR_lin(j,b,i) = Pr_Linear(j,b,i) / ((10^((q(1)-30)/10)) + inter);
                    SNR_lin(j,b,i) = Pr_Linear(j,b,i) / ((10^((q(1)-30)/10)));
                    TP(j,b,i) = 10^-6 * (BW(1))*log2(1+SNR_lin(j,b,i));
                end

                SINR(j,b,i) = 10.*log10(SINR_lin(j,b,i));

                case DataRate(2)

                if R(j,b,i) > R_limit
                    SINR_lin(j,b,i) = 0;
                    TP(j,b,i) = 0;
                else
                    SINR_lin(j,b,i) = Pr_Linear(j,b,i) / ((10^((q(2)-30)/10)) + inter);
                    SNR_lin(j,b,i) = Pr_Linear(j,b,i) / ((10^((q(2)-30)/10)));
                    TP(j,b,i) = 10^-6 * (BW(2))*log2(1+SNR_lin(j,b,i));
                end

                SINR(j,b,i) = 10.*log10(SINR_lin(j,b,i));

                case DataRate(3)

                if R(j,b,i) > R_limit
                    SINR_lin(j,b,i) = 0;
                    TP(j,b,i) = 0;
                else
                    SINR_lin(j,b,i) = Pr_Linear(j,b,i) / ((10^((q(3)-30)/10)) + inter);
                    SNR_lin(j,b,i) = Pr_Linear(j,b,i) / ((10^((q(3)-30)/10)));
                    TP(j,b,i) = 10^-6 * (BW(3))*log2(1+SNR_lin(j,b,i));
                end

                SINR(j,b,i) = 10.*log10(SINR_lin(j,b,i));

            end

            end
    end

% Association and Coverage Analysis for the Current Solution

Assoc = 0;

for j = 1:1:size(SINR, 3)

    SINR_positioning(:,j) = SINR(:,int16(I(:,:,j)),j);
    TP_positioning(:,j) = TP(:,int16(I(:,:,j)),j);
    R_positioning(:,j) = R(:,int16(I(:,:,j)),j);

end

SINR_max = max(SINR_positioning,[],2);
TP_max = max(TP_positioning,[],2);
R_min = min(R_positioning,[],2);

for k = 1:1:size(Users,1)
     if SINR_max(k,1) >= -10 && TP_max(k,1) > 0 && R_min(k,1) <= R_limit
        
         Assoc = Assoc + 1;
    end

end

for b = 1:1:Pos_points

        R_positioning_altered = R_positioning;
        SINR_positioning_altered = SINR_positioning;
        TP_positioning_altered = TP_positioning;

        R_positioning_altered(:,i) = R(:,b,i);
        SINR_positioning_altered(:,i) = SINR(:,b,i);
        TP_positioning_altered(:,i) = TP(:,b,i);
        
        R_min_altered = min(R_positioning_altered,[],2);

        SINR_max_altered = max(SINR_positioning_altered,[],2);
        TP_max_altered = max(TP_positioning_altered,[],2);

    for k = 1:1:size(Users,1)

        if SINR_max_altered(k,1) >= -10 && R_min_altered(k,1) <= R_limit ...
                && TP_max_altered(k,1) > 0
            
            Assoc_Matrix(1,b,i) = Assoc_Matrix(1,b,i) + 1;
    
        end
    end
end

Mean_SINR = mean(SINR,1);

Mean_SINR_lin = mean(SINR_lin,1);

Mean_Coverage = Assoc_Matrix;

if mean(Mean_Coverage(:,:,i),2) == Mean_Coverage(:,1,i)
    Mean_Coverage(:,1,i) = Mean_Coverage(:,1,i) - 1;
end

Mean_TP = mean(TP,1);

Mean_Demand = mean(Demand_Density,1);

Mean_Interference = mean(Interference,1);

Mean_Interference_lin = mean(Interference_lin,1);


%% Apply AHP to decide the best spot

B = transpose([Mean_SINR_lin(:,:,i); Mean_TP(:,:,i); Mean_Demand(:,:,i); Mean_Coverage(:,:,i)]);

names = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9'};

criteria_names = {'SINR';'Throughput';'Demand per Point'; 'Coverage'};

% AHP function

[h, c, a] = ahp([], C, 'beneficial', B, 'names', names, 'criterianames', criteria_names, 'normalization', 'minmax');

Weights(:,i) = table2array(h);

[User(:,:,i), I(:,:,i)] = max(Weights(:,i));

Sol(:,:,i) = [User(:,:,i), I(:,:,i)];

%% Treat metrics to show averages of best positioning

SINR_positioning_final(:,i) = SINR(:,Best(1,2,i),i);
SINR_lin_positioning_final(:,i) = SINR_lin(:,Best(1,2,i),i);

SINR_max_final = max(SINR_positioning_final,[],2);
SINR_lin_max_final = max(SINR_lin_positioning_final,[],2);

Avg_SINR_final = 10.*log10(mean(SINR_lin_max_final));

Demand_Density_final(:,i) = Demand_Density(:,Best(1,2,i),i);
DD_max_final = max(Demand_Density_final, [],2);
Avg_DD_final = mean(DD_max_final);

TP_positioning_final(:,i) = TP(:,Best(1,2,i),i);
TP_max_final = max(TP_positioning_final,[],2);

TP_eMBB = 0;
eMBB_count = 0;

for j = 1:length(TP_max_final)
    if User_Service(j) == 1
        
        TP_eMBB = TP_eMBB + TP_max_final(j);
        eMBB_count = eMBB_count + 1;
    end
end

Avg_TP_final = TP_eMBB/eMBB_count;


end

    if mean(Sol(1,1,:), 3) > mean(Best(1,1,:), 3) && n >= 2
        Best = Sol;
        % disp(strcat('New optimal solution found in iteration: ', num2str(n)));

    else 
        if mean(Sol(1,1,:), 3) <= mean(Best(1,1,:))
            tol = tol + 1;
        end

        if (tol == 2)
            % disp(strcat('No more optimal solutions available. Iterations: ', num2str(n)));
            break

        end
    end

end

AHP_time = toc;

%% Plotting

if plots == true

% Figure 1: drone central, min. space random placement

    % figure
    % hold on
    % title('SINR (dB), Random Min. Spacing')
    % scatter(Users(:,1), Users(:,2), 20, "filled", "MarkerFaceAlpha", 0.85)
    % %scatter(Users(:,1), Users(:,2), 40, 'pentagram', "filled", "MarkerFaceAlpha", 0.90, "LineWidth", 2);
    % for i = 1:1:N_UAV
    %     scatter(Points_Linear(:,1,i), Points_Linear(:,2,i), 60, '+', 'linewidth', 0.8)
    % end
    % viscircles(Center,Radius, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'red');
    % scatter(Center(:,1),Center(:,2), 400, 'pentagram', 'MarkerEdgeColor', 'red', 'linewidth', 1.5)
    % rectangle('Position', [0 0 X Y])
    % axis([-X/10 X+X/10 -Y/10 Y+Y/10])
    % daspect([1 1 1])
    % xlabel('x (m)')
    % ylabel('y (m)')
    % hold off

    % Plot Figure 2:

    figure
    hold on
    title('SINR (dB)');
    scatter(Users(:,1), Users(:,2), 40, SINR_max_final, "filled", "MarkerFaceAlpha", 0.85)
    %scatter(Users(:,1), Users(:,2), 40, SINR_max_final, 'pentagram', "filled", "MarkerFaceAlpha", 0.90, "LineWidth", 2)%, 5, 'MarkerEdgeColor', 'blue', 'linewidth', 0.2)
    colormap("parula")
    colorbar
    for i = 1:1:N_UAV
        scatter(Points_Linear(:,1,i), Points_Linear(:,2,i), 120, '+', 'linewidth', 0.8)
        viscircles([UAV(Best(:,2,i),1,i), UAV(Best(:,2,i),2,i)], R_limit, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'red');
        scatter(UAV(Best(:,2,i),1,i), UAV(Best(:,2,i),2,i), 400, 'pentagram', 'MarkerEdgeColor', 'red', 'linewidth', 1.5)
    end
    rectangle('Position', [0 0 X Y])
    axis([-X/10 X+X/10 -Y/10 Y+Y/10])
    daspect([1 1 1])
    xlabel('x (m)')
    ylabel('y (m)')
    hold off

end

end