clc
clear

tic

dir = '\\maa1.cc.lut.fi\home\z143799\Documents\Projects\uav-ahp-main\';

nmi = 1000;
DataPoints = 30;
Tol = 1.01;
pop_size = 25;

plots = false;

% CaseSelect:

% 1 for SINR prio
% 2 for TP prio
% 3 for users with higher demand prio

% This will give objective functions based on the criteria weights
% given by AHP in each of those 3 options.

CaseSelect = [1 2 3];

%% Simulation 1:

DataRate1 = [100*10^6, 15*10^6, 7*10^6];

Users1 = load(horzcat(dir,'Users1.txt'));

Demand1 = load(horzcat(dir,'Demand1.txt'));

delta_cov1 = 0.20;
N_users1 = length(Users1);
X1 = 2000;
Y1 = 2000;
h1 = 100; %height in meters

alpha1 = 4.88;
beta1 = 0.43;

[Zcov1,Zcap1,Rcov1,Rcap1,Pt1] = numerology(h1,X1,Y1,N_users1,delta_cov1,DataRate1, alpha1, beta1);

AHP_time1 = zeros(length(CaseSelect),DataPoints);
Network_Metrics_AHP1 = zeros(length(CaseSelect),DataPoints,4);

Network_Metrics_CS1 = zeros(length(CaseSelect),DataPoints,4);
CS_time1 = zeros(length(CaseSelect),DataPoints);

Network_Metrics_PSO1 = zeros(length(CaseSelect),DataPoints,4);
PSO_time1 = zeros(length(CaseSelect),DataPoints);

Network_Metrics_NSGA1 = zeros(length(CaseSelect),DataPoints,4);
NSGA_time1 = zeros(length(CaseSelect),DataPoints);

for i = 1:length(CaseSelect)
    for j = 1:DataPoints

    [~, ~, Network_Metrics_AHP1(i,j,:), AHP_time1(i,j)] = drone_positioning_AHP(CaseSelect(i), Pt1, N_users1, X1, ...
    Y1, h1, DataRate1, Users1, Demand1, Zcov1, Zcap1, Rcov1, Rcap1, nmi, alpha1, beta1, plots);

    [~,~, Network_Metrics_PSO1(i,j,:), PSO_time1(i,j)] = PSO(CaseSelect(i), Pt1, N_users1, X1, Y1, h1, DataRate1, ...
        Users1, Demand1, Zcov1, Zcap1, Rcov1, Rcap1, nmi, Tol, pop_size, alpha1, beta1, Network_Metrics_AHP1(i,j,:), plots);
    
    [~, ~, Network_Metrics_CS1(i,j,:), CS_time1(i,j)] = cs_positioning(CaseSelect(i), Pt1, N_users1, X1, ...
        Y1, h1, DataRate1, Users1, Demand1, Zcov1, Zcap1, Rcov1, Rcap1, nmi, Tol, pop_size, alpha1, beta1, Network_Metrics_AHP1(i,j,:), plots);
    
    [~, ~, Network_Metrics_NSGA1(i,j,:), NSGA_time1(i,j)] = nsga2_positioning(CaseSelect(i), Pt1, N_users1, X1, Y1, h1, DataRate1, ...
        Users1, Demand1, Zcov1, Zcap1, Rcov1, Rcap1, nmi, Tol, pop_size, alpha1, beta1, Network_Metrics_AHP1(i,j,:), plots);
    
    [i,j]
    
    end

    Mean_SINR_AHP1(i) = mean(Network_Metrics_AHP1(i,:,1));
    Mean_DD_AHP1(i) = mean(Network_Metrics_AHP1(i,:,2));
    Mean_Assoc1(i) = mean(Network_Metrics_AHP1(i,:,4))/N_users1;
    Mean_AHP_time1(i) = mean(AHP_time1(i,:));

    Mean_SINR_CS1(i) = mean(Network_Metrics_CS1(i,:,1));
    Mean_DD_CS1(i) = mean(Network_Metrics_CS1(i,:,2));
    Mean_Assoc_CS1(i) = mean(Network_Metrics_CS1(i,:,4))/N_users1;
    Mean_CS_time1(i) = mean(CS_time1(i,:));

    Mean_SINR_PSO1(i) = mean(Network_Metrics_PSO1(i,:,1));
    Mean_DD_PSO1(i) = mean(Network_Metrics_PSO1(i,:,2));
    Mean_Assoc_PSO1(i) = mean(Network_Metrics_PSO1(i,:,4))/N_users1;
    Mean_PSO_time1(i) = mean(PSO_time1(i,:));

    Mean_SINR_NSGA1(i) = mean(Network_Metrics_NSGA1(i,:,1));
    Mean_DD_NSGA1(i) = mean(Network_Metrics_NSGA1(i,:,2));
    Mean_Assoc_NSGA1(i) = mean(Network_Metrics_NSGA1(i,:,4))/N_users1;
    Mean_NSGA_time1(i) = mean(NSGA_time1(i,:));

    Err_Assoc1(i) = std(Network_Metrics_AHP1(i,:,4))/N_users1;
    Err_SINR_AHP1(i) = std(Network_Metrics_AHP1(i,:,1));
    Err_DD_AHP1(i) = std(Network_Metrics_AHP1(i,:,2));
    Err_AHP_time1(i) = std(AHP_time1(i,:));

    Err_Assoc_CS1(i) = std(Network_Metrics_CS1(i,:,4))/N_users1;
    Err_SINR_CS1(i) = std(Network_Metrics_CS1(i,:,1)); 
    Err_DD_CS1(i) = std(Network_Metrics_CS1(i,:,2)); 
    Err_CS_time1(i) = std(CS_time1(i,:));

    Err_Assoc_PSO1(i) = std(Network_Metrics_PSO1(i,:,4))/N_users1;
    Err_SINR_PSO1(i) = std(Network_Metrics_PSO1(i,:,1)); 
    Err_DD_PSO1(i) = std(Network_Metrics_PSO1(i,:,2)); 
    Err_PSO_time1(i) = std(PSO_time1(i,:));

    Err_Assoc_NSGA1(i) = std(Network_Metrics_NSGA1(i,:,4))/N_users1;
    Err_SINR_NSGA1(i) = std(Network_Metrics_NSGA1(i,:,1)); 
    Err_DD_NSGA1(i) = std(Network_Metrics_NSGA1(i,:,2)); 
    Err_NSGA_time1(i) = std(NSGA_time1(i,:));
end

%% Simulation 2

DataRate2 = [100*10^6, 15*10^6, 7*10^6];

Users2 = load(horzcat(dir,'Users2.txt'));

Demand2 = load(horzcat(dir,'Demand2.txt'));

delta_cov2 = 0.40;
N_users2 = length(Users2);
X2 = 200;
Y2 = 200;
h2 = 30; %height in meters

alpha2 = 27.23;
beta2 = 0.08;

[Zcov2,Zcap2,Rcov2,Rcap2,Pt2] = numerology(h2,X2,Y2,N_users2,delta_cov2,DataRate2, alpha2, beta2);

AHP_time2 = zeros(length(CaseSelect),DataPoints);
Network_Metrics_AHP2 = zeros(length(CaseSelect),DataPoints,4);

Network_Metrics_CS2 = zeros(length(CaseSelect),DataPoints,4);
CS_time2 = zeros(length(CaseSelect),DataPoints);

Network_Metrics_PSO2 = zeros(length(CaseSelect),DataPoints,4);
PSO_time2 = zeros(length(CaseSelect),DataPoints);

Network_Metrics_NSGA2 = zeros(length(CaseSelect),DataPoints,4);
NSGA_time2 = zeros(length(CaseSelect),DataPoints);

for i = 1:length(CaseSelect)
    for j = 1:DataPoints

    [~, ~, Network_Metrics_AHP2(i,j,:), AHP_time2(i,j)] = drone_positioning_AHP(CaseSelect(i), Pt2, N_users2, X2, ...
    Y2, h2, DataRate2, Users2, Demand2, Zcov2, Zcap2, Rcov2, Rcap2, nmi, alpha2, beta2, plots);

    [~,~, Network_Metrics_PSO2(i,j,:), PSO_time2(i,j)] = PSO(CaseSelect(i), Pt2, N_users2, X2, Y2, h2, DataRate2, ...
        Users2, Demand2, Zcov2, Zcap2, Rcov2, Rcap2, nmi, Tol, pop_size, alpha2, beta2, Network_Metrics_AHP2(i,j,:), plots);
    
    [~, ~, Network_Metrics_CS2(i,j,:), CS_time2(i,j)] = cs_positioning(CaseSelect(i), Pt2, N_users2, X2, ...
        Y2, h2, DataRate2, Users2, Demand2, Zcov2, Zcap2, Rcov2, Rcap2, nmi, Tol, pop_size, alpha2, beta2, Network_Metrics_AHP2(i,j,:), plots);
    
    [~, ~, Network_Metrics_NSGA2(i,j,:), NSGA_time2(i,j)] = nsga2_positioning(CaseSelect(i), Pt2, N_users2, X2, Y2, h2, DataRate2, ...
        Users2, Demand2, Zcov2, Zcap2, Rcov2, Rcap2, nmi, Tol, pop_size, alpha2, beta2, Network_Metrics_AHP2(i,j,:), plots);
    
    [i,j]
    
    end

    Mean_SINR_AHP2(i) = mean(Network_Metrics_AHP2(i,:,1));
    Mean_DD_AHP2(i) = mean(Network_Metrics_AHP2(i,:,2));
    Mean_Assoc2(i) = mean(Network_Metrics_AHP2(i,:,4))/N_users2;
    Mean_AHP_time2(i) = mean(AHP_time2(i,:));

    Mean_SINR_CS2(i) = mean(Network_Metrics_CS2(i,:,1));
    Mean_DD_CS2(i) = mean(Network_Metrics_CS2(i,:,2));
    Mean_Assoc_CS2(i) = mean(Network_Metrics_CS2(i,:,4))/N_users2;
    Mean_CS_time2(i) = mean(CS_time2(i,:));

    Mean_SINR_PSO2(i) = mean(Network_Metrics_PSO2(i,:,1));
    Mean_DD_PSO2(i) = mean(Network_Metrics_PSO2(i,:,2));
    Mean_Assoc_PSO2(i) = mean(Network_Metrics_PSO2(i,:,4))/N_users2;
    Mean_PSO_time2(i) = mean(PSO_time2(i,:));

    Mean_SINR_NSGA2(i) = mean(Network_Metrics_NSGA2(i,:,1));
    Mean_DD_NSGA2(i) = mean(Network_Metrics_NSGA2(i,:,2));
    Mean_Assoc_NSGA2(i) = mean(Network_Metrics_NSGA2(i,:,4))/N_users2;
    Mean_NSGA_time2(i) = mean(NSGA_time2(i,:));

    Err_Assoc2(i) = std(Network_Metrics_AHP2(i,:,4))/N_users2;
    Err_SINR_AHP2(i) = std(Network_Metrics_AHP2(i,:,1));
    Err_DD_AHP2(i) = std(Network_Metrics_AHP2(i,:,2));
    Err_AHP_time2(i) = std(AHP_time2(i,:));

    Err_Assoc_CS2(i) = std(Network_Metrics_CS2(i,:,4))/N_users2;
    Err_SINR_CS2(i) = std(Network_Metrics_CS2(i,:,1)); 
    Err_DD_CS2(i) = std(Network_Metrics_CS2(i,:,2)); 
    Err_CS_time2(i) = std(CS_time2(i,:));

    Err_Assoc_PSO2(i) = std(Network_Metrics_PSO2(i,:,4))/N_users2;
    Err_SINR_PSO2(i) = std(Network_Metrics_PSO2(i,:,1)); 
    Err_DD_PSO2(i) = std(Network_Metrics_PSO2(i,:,2)); 
    Err_PSO_time2(i) = std(PSO_time2(i,:));

    Err_Assoc_NSGA2(i) = std(Network_Metrics_NSGA2(i,:,4))/N_users2;
    Err_SINR_NSGA2(i) = std(Network_Metrics_NSGA2(i,:,1)); 
    Err_DD_NSGA2(i) = std(Network_Metrics_NSGA2(i,:,2)); 
    Err_NSGA_time2(i) = std(NSGA_time2(i,:));
end

%% Simulation 3

DataRate3 = [6.0*10^6 3.0*10^6 0.8*10^6];

Users3 = load(horzcat(dir,'Users_festival.txt'));

Demand3 = load(horzcat(dir,'Demand_festival.txt'));

delta_cov3 = 0.40;
N_users3 = length(Users3);
X3 = 1270;
Y3 = 400;
h3 = 50; %height in meters

alpha3 = 9.6;
beta3 = 0.28;

[Zcov3,Zcap3,Rcov3,Rcap3,Pt3] = numerology(h3,X3,Y3,N_users3,delta_cov3,DataRate3, alpha3, beta3);

AHP_time3 = zeros(length(CaseSelect),DataPoints);
Network_Metrics_AHP3 = zeros(length(CaseSelect),DataPoints,4);

Network_Metrics_CS3 = zeros(length(CaseSelect),DataPoints,4);
CS_time3 = zeros(length(CaseSelect),DataPoints);

Network_Metrics_PSO3 = zeros(length(CaseSelect),DataPoints,4);
PSO_time3 = zeros(length(CaseSelect),DataPoints);

Network_Metrics_NSGA3 = zeros(length(CaseSelect),DataPoints,4);
NSGA_time3 = zeros(length(CaseSelect),DataPoints);

for i = 1:length(CaseSelect)
    for j = 1:DataPoints

    [~, ~, Network_Metrics_AHP3(i,j,:), AHP_time3(i,j)] = drone_positioning_AHP(CaseSelect(i), Pt3, N_users3, X3, ...
    Y3, h3, DataRate3, Users3, Demand3, Zcov3, Zcap3, Rcov3, Rcap3, nmi, alpha3, beta3, plots);

    [~,~, Network_Metrics_PSO3(i,j,:), PSO_time3(i,j)] = PSO(CaseSelect(i), Pt3, N_users3, X3, Y3, h3, DataRate3, ...
        Users3, Demand3, Zcov3, Zcap3, Rcov3, Rcap3, nmi, Tol, pop_size, alpha3, beta3, Network_Metrics_AHP3(i,j,:), plots);
    
    [~, ~, Network_Metrics_CS3(i,j,:), CS_time3(i,j)] = cs_positioning(CaseSelect(i), Pt3, N_users3, X3, ...
        Y3, h3, DataRate3, Users3, Demand3, Zcov3, Zcap3, Rcov3, Rcap3, nmi, Tol, pop_size, alpha3, beta3, Network_Metrics_AHP3(i,j,:), plots);
    
    [~, ~, Network_Metrics_NSGA3(i,j,:), NSGA_time3(i,j)] = nsga2_positioning(CaseSelect(i), Pt3, N_users3, X3, Y3, h3, DataRate3, ...
        Users3, Demand3, Zcov3, Zcap3, Rcov3, Rcap3, nmi, Tol, pop_size, alpha3, beta3, Network_Metrics_AHP3(i,j,:), plots);
    
    [i,j]
    
    end

    Mean_SINR_AHP3(i) = mean(Network_Metrics_AHP3(i,:,1));
    Mean_DD_AHP3(i) = mean(Network_Metrics_AHP3(i,:,2));
    Mean_Assoc3(i) = mean(Network_Metrics_AHP3(i,:,4))/N_users3;
    Mean_AHP_time3(i) = mean(AHP_time3(i,:));

    Mean_SINR_CS3(i) = mean(Network_Metrics_CS3(i,:,1));
    Mean_DD_CS3(i) = mean(Network_Metrics_CS3(i,:,2));
    Mean_Assoc_CS3(i) = mean(Network_Metrics_CS3(i,:,4))/N_users3;
    Mean_CS_time3(i) = mean(CS_time3(i,:));

    Mean_SINR_PSO3(i) = mean(Network_Metrics_PSO3(i,:,1));
    Mean_DD_PSO3(i) = mean(Network_Metrics_PSO3(i,:,2));
    Mean_Assoc_PSO3(i) = mean(Network_Metrics_PSO3(i,:,4))/N_users3;
    Mean_PSO_time3(i) = mean(PSO_time3(i,:));

    Mean_SINR_NSGA3(i) = mean(Network_Metrics_NSGA3(i,:,1));
    Mean_DD_NSGA3(i) = mean(Network_Metrics_NSGA3(i,:,2));
    Mean_Assoc_NSGA3(i) = mean(Network_Metrics_NSGA3(i,:,4))/N_users3;
    Mean_NSGA_time3(i) = mean(NSGA_time3(i,:));

    Err_Assoc3(i) = std(Network_Metrics_AHP3(i,:,4))/N_users3;
    Err_SINR_AHP3(i) = std(Network_Metrics_AHP3(i,:,1));
    Err_DD_AHP3(i) = std(Network_Metrics_AHP3(i,:,2));
    Err_AHP_time3(i) = std(AHP_time3(i,:));

    Err_Assoc_CS3(i) = std(Network_Metrics_CS3(i,:,4))/N_users3;
    Err_SINR_CS3(i) = std(Network_Metrics_CS3(i,:,1)); 
    Err_DD_CS3(i) = std(Network_Metrics_CS3(i,:,2)); 
    Err_CS_time3(i) = std(CS_time3(i,:));

    Err_Assoc_PSO3(i) = std(Network_Metrics_PSO3(i,:,4))/N_users3;
    Err_SINR_PSO3(i) = std(Network_Metrics_PSO3(i,:,1)); 
    Err_DD_PSO3(i) = std(Network_Metrics_PSO3(i,:,2)); 
    Err_PSO_time3(i) = std(PSO_time3(i,:));

    Err_Assoc_NSGA3(i) = std(Network_Metrics_NSGA3(i,:,4))/N_users3;
    Err_SINR_NSGA3(i) = std(Network_Metrics_NSGA3(i,:,1)); 
    Err_DD_NSGA3(i) = std(Network_Metrics_NSGA3(i,:,2)); 
    Err_NSGA_time3(i) = std(NSGA_time3(i,:));
end



%% Plots

ErrData1 = [Err_Assoc_NSGA3(3), Err_Assoc_CS3(3), Err_Assoc_PSO3(3), Err_Assoc3(3);
    Err_Assoc_NSGA3(2), Err_Assoc_CS3(2), Err_Assoc_PSO3(2), Err_Assoc3(2);
    Err_Assoc_NSGA3(1), Err_Assoc_CS3(1), Err_Assoc_PSO3(1), Err_Assoc3(1);
    Err_Assoc_NSGA2(3), Err_Assoc_CS2(3), Err_Assoc_PSO2(3), Err_Assoc2(3);
    Err_Assoc_NSGA2(2), Err_Assoc_CS2(2), Err_Assoc_PSO2(2), Err_Assoc2(2);
    Err_Assoc_NSGA2(1), Err_Assoc_CS2(1), Err_Assoc_PSO2(1), Err_Assoc2(1);
    Err_Assoc_NSGA1(3), Err_Assoc_CS1(3), Err_Assoc_PSO1(3), Err_Assoc1(3);
    Err_Assoc_NSGA1(2), Err_Assoc_CS1(2), Err_Assoc_PSO1(2), Err_Assoc1(2);
    Err_Assoc_NSGA1(1), Err_Assoc_CS1(1), Err_Assoc_PSO1(1), Err_Assoc1(1)
];

PlotData1 = [Mean_Assoc_NSGA3(3), Mean_Assoc_CS3(3), Mean_Assoc_PSO3(3), Mean_Assoc3(3);
    Mean_Assoc_NSGA3(2), Mean_Assoc_CS3(2), Mean_Assoc_PSO3(2), Mean_Assoc3(2);
    Mean_Assoc_NSGA3(1), Mean_Assoc_CS3(1), Mean_Assoc_PSO3(1), Mean_Assoc3(1);
    Mean_Assoc_NSGA2(3), Mean_Assoc_CS2(3), Mean_Assoc_PSO2(3), Mean_Assoc2(3);
    Mean_Assoc_NSGA2(2), Mean_Assoc_CS2(2), Mean_Assoc_PSO2(2), Mean_Assoc2(2);
    Mean_Assoc_NSGA2(1), Mean_Assoc_CS2(1), Mean_Assoc_PSO2(1), Mean_Assoc2(1);
    Mean_Assoc_NSGA1(3), Mean_Assoc_CS1(3), Mean_Assoc_PSO1(3), Mean_Assoc1(3);
    Mean_Assoc_NSGA1(2), Mean_Assoc_CS1(2), Mean_Assoc_PSO1(2), Mean_Assoc1(2);
    Mean_Assoc_NSGA1(1), Mean_Assoc_CS1(1), Mean_Assoc_PSO1(1), Mean_Assoc1(1)
];

% Plot 1: Coverage

figure('Units', 'pixels', 'Position', [0, 0, 1400, 2000]);
hold on

% === Label Setup ===
numScenarios = 3;
numPlans = 3;
numMethods = 4;

scenarioLabels = {'Scenario 1', 'Scenario 2', 'Scenario 3'};
planLabels = {'SINR', 'TP', 'Coverage'};  % Reverse order from your matrix layout
methodLabels = {'UAV-AHP','PSO', 'CS', 'NSGA-II'};
methodColors = [0, 114, 178; 240, 228, 66; 213, 94, 0; 68, 170, 153]/255;

% === Reverse rows to match order: Scenario 1 (bottom) to 3 (top)
PlotData1 = fliplr(flipud(PlotData1));  % Now SINR1 → Cov3 order
ErrData1 = fliplr(flipud(ErrData1));

% === Layout Configuration ===
methodSpacing = 1;
planSpacing = 5;
scenarioSpacing = 2;

% === Generate Y positions ===
y = [];
yGroupCenters = zeros(numScenarios, numPlans);
barIdx = 0;

for s = 1:numScenarios
    for p = 1:numPlans
        baseY = (s-1)*(numPlans*planSpacing + scenarioSpacing) + (p-1)*planSpacing;
        for m = 1:numMethods
            barIdx = barIdx + 1;
            y(barIdx) = baseY + m * methodSpacing;
        end
        yGroupCenters(s, p) = baseY + (numMethods + 1)/2 * methodSpacing;
    end
end

% === Flatten data for plotting ===
flatData = reshape(PlotData1', [], 1);
flatErrors = reshape(ErrData1', [], 1);

% === Plot Bars and Errors ===
barHandles = nan(numMethods, 1);  % Preallocate for legend handles

for i = 1:length(flatData)
    methodIdx = mod(i-1, numMethods) + 1;

    % Plot bar and save handle if first time for this method
    hBar = barh(y(i), flatData(i), 0.9, ...
        'FaceColor', methodColors(methodIdx,:), 'EdgeColor', 'k', 'LineWidth', 0.05);

    if isnan(barHandles(methodIdx))
        barHandles(methodIdx) = hBar;
    end

    % Plot error bar without adding to legend
    errorbar(flatData(i), y(i), flatErrors(i), 'horizontal', ...
        'Color', 'k', 'LineWidth', 0.2, 'HandleVisibility', 'off');
end

% === Dashed Lines between Scenarios ===
maxY = max(y) + 3;
for s = 1:numScenarios-1
    splitY = s * (numPlans * planSpacing + scenarioSpacing) - scenarioSpacing / 2;
    plot(xlim, [splitY splitY], 'k--', 'LineWidth', 0.6);
end

% === Y-Axis Plan Labels ===
yticks(reshape(yGroupCenters', [], 1));
yticklabels(repmat(planLabels, 1, numScenarios));

% === X-Axis Settings ===
xlabel('Mean Coverage of Users')
set(gca,'XTickLabel',{'65%','70%','75%','80%','85%','90%','95%','100%'})
xlim([0.65 1.02])
ylim([0 maxY])
set(gca, 'YDir', 'reverse')

% === Add Scenario Labels ===
for s = 1:numScenarios
    % Compute vertical center of this scenario's full group
    groupTop = (s-1)*(numPlans*planSpacing + scenarioSpacing);
    groupBottom = groupTop + numPlans*planSpacing;
    groupCenterY = (groupTop + groupBottom) / 2;

    text(0.635, groupCenterY, scenarioLabels{s}, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 10)
end

% === Legend and Title ===
legend(barHandles, methodLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Coverage Performance by Scenario, Plan, and Method')
box on
hold off

% Plot 2: Running time

ErrData2 = [Err_NSGA_time3(3), Err_CS_time3(3), Err_PSO_time3(3), Err_AHP_time3(3);
    Err_NSGA_time3(2), Err_CS_time3(2), Err_PSO_time3(2), Err_AHP_time3(2);
    Err_NSGA_time3(1), Err_CS_time3(1), Err_PSO_time3(1), Err_AHP_time3(1);
    Err_NSGA_time2(3), Err_CS_time2(3), Err_PSO_time2(3), Err_AHP_time2(3);
    Err_NSGA_time2(2), Err_CS_time2(2), Err_PSO_time2(2), Err_AHP_time2(2);
    Err_NSGA_time2(1), Err_CS_time2(1), Err_PSO_time2(1), Err_AHP_time2(1);
    Err_NSGA_time1(3), Err_CS_time1(3), Err_PSO_time1(3), Err_AHP_time1(3);
    Err_NSGA_time1(2), Err_CS_time1(2), Err_PSO_time1(2), Err_AHP_time1(2);
    Err_NSGA_time1(1), Err_CS_time1(1), Err_PSO_time1(1), Err_AHP_time1(1)
];

PlotData2 = [Mean_NSGA_time3(3), Mean_CS_time3(3), Mean_PSO_time3(3), Mean_AHP_time3(3);
    Mean_NSGA_time3(2), Mean_CS_time3(2), Mean_PSO_time3(2), Mean_AHP_time3(2);
    Mean_NSGA_time3(1), Mean_CS_time3(1), Mean_PSO_time3(1), Mean_AHP_time3(1);
    Mean_NSGA_time2(3), Mean_CS_time2(3), Mean_PSO_time2(3), Mean_AHP_time2(3);
    Mean_NSGA_time2(2), Mean_CS_time2(2), Mean_PSO_time2(2), Mean_AHP_time2(2);
    Mean_NSGA_time2(1), Mean_CS_time2(1), Mean_PSO_time2(1), Mean_AHP_time2(1);
    Mean_NSGA_time1(3), Mean_CS_time1(3), Mean_PSO_time1(3), Mean_AHP_time1(3);
    Mean_NSGA_time1(2), Mean_CS_time1(2), Mean_PSO_time1(2), Mean_AHP_time1(2);
    Mean_NSGA_time1(1), Mean_CS_time1(1), Mean_PSO_time1(1), Mean_AHP_time1(1)
];

figure('Units', 'pixels', 'Position', [0, 0, 1400, 2000]);
hold on

% === Label Setup ===
numScenarios = 3;
numPlans = 3;
numMethods = 4;

scenarioLabels = {'Scenario 1', 'Scenario 2', 'Scenario 3'};
planLabels = {'SINR', 'TP', 'Coverage'};  % Reverse order from your matrix layout
methodLabels = {'UAV-AHP','PSO', 'CS', 'NSGA-II'};
methodColors = [0, 114, 178; 240, 228, 66; 213, 94, 0; 68, 170, 153]/255;

% === Reverse rows to match order: Scenario 1 (bottom) to 3 (top)
PlotData2 = fliplr(flipud(PlotData2));  % Now SINR1 → Cov3 order
ErrData2 = fliplr(flipud(ErrData2));

% === Layout Configuration ===
methodSpacing = 1;
planSpacing = 5;
scenarioSpacing = 2;

% === Generate Y positions ===
y = [];
yGroupCenters = zeros(numScenarios, numPlans);
barIdx = 0;

for s = 1:numScenarios
    for p = 1:numPlans
        baseY = (s-1)*(numPlans*planSpacing + scenarioSpacing) + (p-1)*planSpacing;
        for m = 1:numMethods
            barIdx = barIdx + 1;
            y(barIdx) = baseY + m * methodSpacing;
        end
        yGroupCenters(s, p) = baseY + (numMethods + 1)/2 * methodSpacing;
    end
end

% === Flatten data for plotting ===
flatData = reshape(PlotData2', [], 1);
flatErrors = reshape(ErrData2', [], 1);

% === Plot Bars and Errors ===
barHandles = nan(numMethods, 1);  % Preallocate for legend handles

for i = 1:length(flatData)
    methodIdx = mod(i-1, numMethods) + 1;

    % Plot bar and save handle if first time for this method
    hBar = barh(y(i), flatData(i), 0.9, ...
        'FaceColor', methodColors(methodIdx,:), 'EdgeColor', 'k', 'LineWidth', 0.05);

    if isnan(barHandles(methodIdx))
        barHandles(methodIdx) = hBar;
    end

    % Plot error bar without adding to legend
    errorbar(flatData(i), y(i), flatErrors(i), 'horizontal', ...
        'Color', 'k', 'LineWidth', 0.2, 'HandleVisibility', 'off');
end

% === Dashed Lines between Scenarios ===
maxY = max(y) + 3;
for s = 1:numScenarios-1
    splitY = s * (numPlans * planSpacing + scenarioSpacing) - scenarioSpacing / 2;
    plot(xlim, [splitY splitY], 'k--', 'LineWidth', 0.6);
end

% === Y-Axis Plan Labels ===
yticks(reshape(yGroupCenters', [], 1));
yticklabels(repmat(planLabels, 1, numScenarios));

% === X-Axis Settings ===
xlabel('Mean Running Time (seconds)')
% set(gca,'XTickLabel',{'65%','70%','75%','80%','85%','90%','95%','100%'})
% xlim([0.65 1.02])
xscale("log")
ylim([0 maxY])
set(gca, 'YDir', 'reverse')

% === Add Scenario Labels ===
for s = 1:numScenarios
    % Compute vertical center of this scenario's full group
    groupTop = (s-1)*(numPlans*planSpacing + scenarioSpacing);
    groupBottom = groupTop + numPlans*planSpacing;
    groupCenterY = (groupTop + groupBottom) / 2;

    text(0.006, groupCenterY, scenarioLabels{s}, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 10)
end

% === Legend and Title ===
legend(barHandles, methodLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Running Time by Scenario, Plan, and Method')
box on
hold off

% Plot 3: SINR

ErrData3 = [Err_SINR_NSGA3(3), Err_SINR_CS3(3), Err_SINR_PSO3(3), Err_SINR_AHP3(3);
    Err_SINR_NSGA3(2), Err_SINR_CS3(2), Err_SINR_PSO3(2), Err_SINR_AHP3(2);
    Err_SINR_NSGA3(1), Err_SINR_CS3(1), Err_SINR_PSO3(1), Err_SINR_AHP3(1);
    Err_SINR_NSGA2(3), Err_SINR_CS2(3), Err_SINR_PSO2(3), Err_SINR_AHP2(3);
    Err_SINR_NSGA2(2), Err_SINR_CS2(2), Err_SINR_PSO2(2), Err_SINR_AHP2(2);
    Err_SINR_NSGA2(1), Err_SINR_CS2(1), Err_SINR_PSO2(1), Err_SINR_AHP2(1);
    Err_SINR_NSGA1(3), Err_SINR_CS1(3), Err_SINR_PSO1(3), Err_SINR_AHP1(3);
    Err_SINR_NSGA1(2), Err_SINR_CS1(2), Err_SINR_PSO1(2), Err_SINR_AHP1(2);
    Err_SINR_NSGA1(1), Err_SINR_CS1(1), Err_SINR_PSO1(1), Err_SINR_AHP1(1);
];

PlotData3 = [Mean_SINR_NSGA3(3), Mean_SINR_CS3(3), Mean_SINR_PSO3(3), Mean_SINR_AHP3(3);
    Mean_SINR_NSGA3(2), Mean_SINR_CS3(2), Mean_SINR_PSO3(2), Mean_SINR_AHP3(2);
    Mean_SINR_NSGA3(1), Mean_SINR_CS3(1), Mean_SINR_PSO3(1), Mean_SINR_AHP3(1);
    Mean_SINR_NSGA2(3), Mean_SINR_CS2(3), Mean_SINR_PSO2(3), Mean_SINR_AHP2(3);
    Mean_SINR_NSGA2(2), Mean_SINR_CS2(2), Mean_SINR_PSO2(2), Mean_SINR_AHP2(2);
    Mean_SINR_NSGA2(1), Mean_SINR_CS2(1), Mean_SINR_PSO2(1), Mean_SINR_AHP2(1);
    Mean_SINR_NSGA1(3), Mean_SINR_CS1(3), Mean_SINR_PSO1(3), Mean_SINR_AHP1(3);
    Mean_SINR_NSGA1(2), Mean_SINR_CS1(2), Mean_SINR_PSO1(2), Mean_SINR_AHP1(2);
    Mean_SINR_NSGA1(1), Mean_SINR_CS1(1), Mean_SINR_PSO1(1), Mean_SINR_AHP1(1);
];

figure('Units', 'pixels', 'Position', [0, 0, 1400, 2000]);
hold on

% === Label Setup ===
numScenarios = 3;
numPlans = 3;
numMethods = 4;

scenarioLabels = {'Scenario 1', 'Scenario 2', 'Scenario 3'};
planLabels = {'SINR', 'TP', 'Coverage'};  % Reverse order from your matrix layout
methodLabels = {'UAV-AHP','PSO', 'CS', 'NSGA-II'};
methodColors = [0, 114, 178; 240, 228, 66; 213, 94, 0; 68, 170, 153]/255;

% === Reverse rows to match order: Scenario 1 (bottom) to 3 (top)
PlotData3 = fliplr(flipud(PlotData3));  % Now SINR1 → Cov3 order
ErrData3 = fliplr(flipud(ErrData3));

% === Layout Configuration ===
methodSpacing = 1;
planSpacing = 5;
scenarioSpacing = 2;

% === Generate Y positions ===
y = [];
yGroupCenters = zeros(numScenarios, numPlans);
barIdx = 0;

for s = 1:numScenarios
    for p = 1:numPlans
        baseY = (s-1)*(numPlans*planSpacing + scenarioSpacing) + (p-1)*planSpacing;
        for m = 1:numMethods
            barIdx = barIdx + 1;
            y(barIdx) = baseY + m * methodSpacing;
        end
        yGroupCenters(s, p) = baseY + (numMethods + 1)/2 * methodSpacing;
    end
end

% === Flatten data for plotting ===
flatData = reshape(PlotData3', [], 1);
flatErrors = reshape(ErrData3', [], 1);

% === Plot Bars and Errors ===
barHandles = nan(numMethods, 1);  % Preallocate for legend handles

for i = 1:length(flatData)
    methodIdx = mod(i-1, numMethods) + 1;

    % Plot bar and save handle if first time for this method
    hBar = barh(y(i), flatData(i), 0.9, ...
        'FaceColor', methodColors(methodIdx,:), 'EdgeColor', 'k', 'LineWidth', 0.05);

    if isnan(barHandles(methodIdx))
        barHandles(methodIdx) = hBar;
    end

    % Plot error bar without adding to legend
    errorbar(flatData(i), y(i), flatErrors(i), 'horizontal', ...
        'Color', 'k', 'LineWidth', 0.2, 'HandleVisibility', 'off');
end

% === Dashed Lines between Scenarios ===
maxY = max(y) + 3;
for s = 1:numScenarios-1
    splitY = s * (numPlans * planSpacing + scenarioSpacing) - scenarioSpacing / 2;
    plot(xlim, [splitY splitY], 'k--', 'LineWidth', 0.6);
end

% === Y-Axis Plan Labels ===
yticks(reshape(yGroupCenters', [], 1));
yticklabels(repmat(planLabels, 1, numScenarios));

% === X-Axis Settings ===
xlabel('Mean SINR (dB)')
% set(gca,'XTickLabel',{'65%','70%','75%','80%','85%','90%','95%','100%'})
% xlim([0.65 1.02])
ylim([0 maxY])
set(gca, 'YDir', 'reverse')

% === Add Scenario Labels ===
for s = 1:numScenarios
    % Compute vertical center of this scenario's full group
    groupTop = (s-1)*(numPlans*planSpacing + scenarioSpacing);
    groupBottom = groupTop + numPlans*planSpacing;
    groupCenterY = (groupTop + groupBottom) / 2;

    text(-1, groupCenterY, scenarioLabels{s}, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 10)
end

% === Legend and Title ===
legend(barHandles, methodLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('SINR Performance by Scenario, Plan, and Method')
box on
hold off

% Plot 4: TP in eMBB slice

ErrData4 = [Err_DD_NSGA3(3), Err_DD_CS3(3), Err_DD_PSO3(3), Err_DD_AHP3(3);
    Err_DD_NSGA3(2), Err_DD_CS3(2), Err_DD_PSO3(2), Err_DD_AHP3(2);
    Err_DD_NSGA3(1), Err_DD_CS3(1), Err_DD_PSO3(1), Err_DD_AHP3(1);
    Err_DD_NSGA2(3), Err_DD_CS2(3), Err_DD_PSO2(3), Err_DD_AHP2(3);
    Err_DD_NSGA2(2), Err_DD_CS2(2), Err_DD_PSO2(2), Err_DD_AHP2(2)
    Err_DD_NSGA2(1), Err_DD_CS2(1), Err_DD_PSO2(1), Err_DD_AHP2(1);
    Err_DD_NSGA1(3), Err_DD_CS1(3), Err_DD_PSO1(3), Err_DD_AHP1(3);
    Err_DD_NSGA1(2), Err_DD_CS1(2), Err_DD_PSO1(2), Err_DD_AHP1(2);
    Err_DD_NSGA1(1), Err_DD_CS1(1), Err_DD_PSO1(1), Err_DD_AHP1(1);
];

PlotData4 = [Mean_DD_NSGA3(3), Mean_DD_CS3(3), Mean_DD_PSO3(3), Mean_DD_AHP3(3);
    Mean_DD_NSGA3(2), Mean_DD_CS3(2), Mean_DD_PSO3(2), Mean_DD_AHP3(2);
    Mean_DD_NSGA3(1), Mean_DD_CS3(1), Mean_DD_PSO3(1), Mean_DD_AHP3(1);
    Mean_DD_NSGA2(3), Mean_DD_CS2(3), Mean_DD_PSO2(3), Mean_DD_AHP2(3);
    Mean_DD_NSGA2(2), Mean_DD_CS2(2), Mean_DD_PSO2(2), Mean_DD_AHP2(2)
    Mean_DD_NSGA2(1), Mean_DD_CS2(1), Mean_DD_PSO2(1), Mean_DD_AHP2(1);
    Mean_DD_NSGA1(3), Mean_DD_CS1(3), Mean_DD_PSO1(3), Mean_DD_AHP1(3);
    Mean_DD_NSGA1(2), Mean_DD_CS1(2), Mean_DD_PSO1(2), Mean_DD_AHP1(2);
    Mean_DD_NSGA1(1), Mean_DD_CS1(1), Mean_DD_PSO1(1), Mean_DD_AHP1(1);
];

figure('Units', 'pixels', 'Position', [0, 0, 1400, 2000]);
hold on

% === Label Setup ===
numScenarios = 3;
numPlans = 3;
numMethods = 4;

scenarioLabels = {'Scenario 1', 'Scenario 2', 'Scenario 3'};
planLabels = {'SINR', 'TP', 'Coverage'};  % Reverse order from your matrix layout
methodLabels = {'UAV-AHP','PSO', 'CS', 'NSGA-II'};
methodColors = [0, 114, 178; 240, 228, 66; 213, 94, 0; 68, 170, 153]/255;

% === Reverse rows to match order: Scenario 1 (bottom) to 3 (top)
PlotData4 = fliplr(flipud(PlotData4));  % Now SINR1 → Cov3 order
ErrData4 = fliplr(flipud(ErrData4));

% === Layout Configuration ===
methodSpacing = 1;
planSpacing = 5;
scenarioSpacing = 2;

% === Generate Y positions ===
y = [];
yGroupCenters = zeros(numScenarios, numPlans);
barIdx = 0;

for s = 1:numScenarios
    for p = 1:numPlans
        baseY = (s-1)*(numPlans*planSpacing + scenarioSpacing) + (p-1)*planSpacing;
        for m = 1:numMethods
            barIdx = barIdx + 1;
            y(barIdx) = baseY + m * methodSpacing;
        end
        yGroupCenters(s, p) = baseY + (numMethods + 1)/2 * methodSpacing;
    end
end

% === Flatten data for plotting ===
flatData = reshape(PlotData4', [], 1);
flatErrors = reshape(ErrData4', [], 1);

% === Plot Bars and Errors ===
barHandles = nan(numMethods, 1);  % Preallocate for legend handles

for i = 1:length(flatData)
    methodIdx = mod(i-1, numMethods) + 1;

    % Plot bar and save handle if first time for this method
    hBar = barh(y(i), flatData(i), 0.9, ...
        'FaceColor', methodColors(methodIdx,:), 'EdgeColor', 'k', 'LineWidth', 0.05);

    if isnan(barHandles(methodIdx))
        barHandles(methodIdx) = hBar;
    end

    % Plot error bar without adding to legend
    errorbar(flatData(i), y(i), flatErrors(i), 'horizontal', ...
        'Color', 'k', 'LineWidth', 0.2, 'HandleVisibility', 'off');
end

% === Dashed Lines between Scenarios ===
maxY = max(y) + 3;
for s = 1:numScenarios-1
    splitY = s * (numPlans * planSpacing + scenarioSpacing) - scenarioSpacing / 2;
    plot(xlim, [splitY splitY], 'k--', 'LineWidth', 0.6);
end

% === Y-Axis Plan Labels ===
yticks(reshape(yGroupCenters', [], 1));
yticklabels(repmat(planLabels, 1, numScenarios));

% === X-Axis Settings ===
xlabel('Mean Throughput in eMBB Slice (Mbits/second)')
% set(gca,'XTickLabel',{'65%','70%','75%','80%','85%','90%','95%','100%'})
xlim([300 700])
ylim([0 maxY])
set(gca, 'YDir', 'reverse')

% === Add Scenario Labels ===
for s = 1:numScenarios
    % Compute vertical center of this scenario's full group
    groupTop = (s-1)*(numPlans*planSpacing + scenarioSpacing);
    groupBottom = groupTop + numPlans*planSpacing;
    groupCenterY = (groupTop + groupBottom) / 2;

    text(285, groupCenterY, scenarioLabels{s}, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 10)
end

% === Legend and Title ===
legend(barHandles, methodLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Throughput Performance by Scenario, Plan, and Method')
box on
hold off

Plot_time = toc;