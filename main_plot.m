clc
clear

tic

dir = '/home/henry-ferreira/Documentos/WORK/PeerJ_2/Trab final/';

nmi = 1000;
DataPoints = 30;

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

alpha = 4.88;
beta = 0.43;

Tol = 1;

[Zcov1,Zcap1,Rcov1,Rcap1,Pt1] = numerology(h1,X1,Y1,N_users1,delta_cov1,DataRate1, alpha, beta);

AHP_time1 = zeros(length(CaseSelect),DataPoints);
Assoc1 = zeros(length(CaseSelect),DataPoints);
Avg_SINR_AHP1 = zeros(length(CaseSelect),DataPoints);
Avg_DD_AHP1 = zeros(length(CaseSelect),DataPoints);

fmin1 = zeros(length(CaseSelect),DataPoints);
Assoc_CS1 = zeros(length(CaseSelect),DataPoints);
Avg_SINR_CS1 = zeros(length(CaseSelect),DataPoints);
Avg_DD_CS1 = zeros(length(CaseSelect),DataPoints);
CS_time1 = zeros(length(CaseSelect),DataPoints);

for i = 1:length(CaseSelect)
    for j = 1:DataPoints
    
    [~, ~, AHP_time1(i,j), Assoc1(i,j), Avg_SINR_AHP1(i,j), Avg_DD_AHP1(i,j)] = drone_positioning_AHP(CaseSelect(i),Pt1, N_users1, X1, ...
        Y1, h1, DataRate1, Users1, Demand1, Zcov1, Zcap1, Rcov1, Rcap1, nmi, alpha, beta, plots);
    
    [~, fmin1(i,j), Assoc_CS1(i,j), Avg_SINR_CS1(i,j), Avg_DD_CS1(i,j), CS_time1(i,j)]=cs_positioning(CaseSelect(i),Pt1, N_users1, X1, ...
        Y1, h1, DataRate1, Users1, Demand1, Zcov1, Zcap1, Rcov1, Rcap1, nmi, Tol, alpha, beta, plots);
    
    end

    Mean_Assoc1(i) = mean(Assoc1(i,:))/N_users1;
    Mean_SINR_AHP1(i) = mean(Avg_SINR_AHP1(i,:));
    Mean_DD_AHP1(i) = mean(Avg_DD_AHP1(i,:));
    Mean_AHP_time1(i) = mean(AHP_time1(i,:));
    
    Mean_Assoc_CS1(i) = mean(Assoc_CS1(i,:))/N_users1;
    Mean_SINR_CS1(i) = mean(Avg_SINR_CS1(i,:)); 
    Mean_DD_CS1(i) = mean(Avg_DD_CS1(i,:)); 
    Mean_CS_time1(i) = mean(CS_time1(i,:));

    Err_Assoc1(i) = std(Assoc1(i,:))/N_users1;
    Err_SINR_AHP1(i) = std(Avg_SINR_AHP1(i,:));
    Err_DD_AHP1(i) = std(Avg_DD_AHP1(i,:));
    Err_AHP_time1(i) = std(AHP_time1(i,:));

    Err_Assoc_CS1(i) = std(Assoc_CS1(i,:))/N_users1;
    Err_SINR_CS1(i) = std(Avg_SINR_CS1(i,:)); 
    Err_DD_CS1(i) = std(Avg_DD_CS1(i,:)); 
    Err_CS_time1(i) = std(CS_time1(i,:));
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

alpha = 27.23;
beta = 0.08;

Tol = 1;

[Zcov2,Zcap2,Rcov2,Rcap2,Pt2] = numerology(h2,X2,Y2,N_users2,delta_cov2,DataRate2, alpha, beta);

AHP_time2 = zeros(length(CaseSelect),DataPoints);
Assoc2 = zeros(length(CaseSelect),DataPoints);
Avg_SINR_AHP2 = zeros(length(CaseSelect),DataPoints);
Avg_DD_AHP2 = zeros(length(CaseSelect),DataPoints);

fmin2 = zeros(length(CaseSelect),DataPoints);
Assoc_CS2 = zeros(length(CaseSelect),DataPoints);
Avg_SINR_CS2 = zeros(length(CaseSelect),DataPoints);
Avg_DD_CS2 = zeros(length(CaseSelect),DataPoints);
CS_time2 = zeros(length(CaseSelect),DataPoints);

for i = 1:length(CaseSelect)
    for j = 1:DataPoints
    
    [~, ~, AHP_time2(i,j), Assoc2(i,j), Avg_SINR_AHP2(i,j), Avg_DD_AHP2(i,j)] = drone_positioning_AHP(CaseSelect(i),Pt2, N_users2, X2, ...
        Y2, h2, DataRate2, Users2, Demand2, Zcov2, Zcap2, Rcov2, Rcap2, nmi, alpha, beta, plots);
    
    [~, fmin2(i,j), Assoc_CS2(i,j), Avg_SINR_CS2(i,j), Avg_DD_CS2(i,j), CS_time2(i,j)]=cs_positioning(CaseSelect(i),Pt2, N_users2, X2, ...
        Y2, h2, DataRate2, Users2, Demand2, Zcov2, Zcap2, Rcov2, Rcap2, nmi, Tol, alpha, beta, plots);
    
    end

    Mean_Assoc2(i) = mean(Assoc2(i,:))/N_users2;
    Mean_SINR_AHP2(i) = mean(Avg_SINR_AHP2(i,:));
    Mean_DD_AHP2(i) = mean(Avg_DD_AHP2(i,:));
    Mean_AHP_time2(i) = mean(AHP_time2(i,:));
    
    Mean_Assoc_CS2(i) = mean(Assoc_CS2(i,:))/N_users2;
    Mean_SINR_CS2(i) = mean(Avg_SINR_CS2(i,:));
    Mean_DD_CS2(i) = mean(Avg_DD_CS2(i,:));
    Mean_CS_time2(i) = mean(CS_time2(i,:));

    Err_Assoc2(i) = std(Assoc2(i,:))/N_users2;
    Err_SINR_AHP2(i) = std(Avg_SINR_AHP2(i,:));
    Err_DD_AHP2(i) = std(Avg_DD_AHP2(i,:));
    Err_AHP_time2(i) = std(AHP_time2(i,:));

    Err_Assoc_CS2(i) = std(Assoc_CS2(i,:))/N_users2;
    Err_SINR_CS2(i) = std(Avg_SINR_CS2(i,:)); 
    Err_DD_CS2(i) = std(Avg_DD_CS2(i,:)); 
    Err_CS_time2(i) = std(CS_time2(i,:));
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

alpha = 9.6;
beta = 0.28;

[Zcov3,Zcap3,Rcov3,Rcap3,Pt3] = numerology(h3,X3,Y3,N_users3,delta_cov3,DataRate3, alpha, beta);

AHP_time3 = zeros(length(CaseSelect),DataPoints);
Assoc3 = zeros(length(CaseSelect),DataPoints);
Avg_SINR_AHP3 = zeros(length(CaseSelect),DataPoints);
Avg_DD_AHP3 = zeros(length(CaseSelect),DataPoints);

fmin3 = zeros(length(CaseSelect),1);
Assoc_CS3 = zeros(length(CaseSelect),1);
Avg_SINR_CS3 = zeros(length(CaseSelect),1);
Avg_DD_CS3 = zeros(length(CaseSelect),1);
CS_time3 = zeros(length(CaseSelect),1);

for i = 1:length(CaseSelect)
    for j = 1:DataPoints

    [~, ~, AHP_time3(i,j), Assoc3(i,j), Avg_SINR_AHP3(i,j), Avg_DD_AHP3(i,j)] = drone_positioning_AHP(CaseSelect(i),Pt3,N_users3, X3, ...
        Y3, h3, DataRate3, Users3, Demand3, Zcov3, Zcap3, Rcov3, Rcap3, nmi, alpha, beta, plots);

    end

    Mean_Assoc3(i) = mean(Assoc3(i,:))/N_users3;
    Mean_SINR_AHP3(i) = mean(Avg_SINR_AHP3(i,:));
    Mean_DD_AHP3(i) = mean(Avg_DD_AHP3(i,:));
    Mean_AHP_time3(i) = mean(AHP_time3(i,:));

    Err_Assoc3(i) = std(Assoc3(i,:))/N_users3;
    Err_SINR_AHP3(i) = std(Avg_SINR_AHP3(i,:));
    Err_DD_AHP3(i) = std(Avg_DD_AHP3(i,:));
    Err_AHP_time3(i) = std(AHP_time3(i,:));

end

parfor i = 1:length(CaseSelect)
    [~, fmin3(i), Assoc_CS3(i), Avg_SINR_CS3(i), Avg_DD_CS3(i), CS_time3(i)]=cs_positioning(CaseSelect(i),Pt3,N_users3, X3, ...
Y3, h3, DataRate3, Users3, Demand3, Zcov3, Zcap3, Rcov3, Rcap3, nmi, Tol, alpha, beta, plots);

    Mean_Assoc_CS3(i) = Assoc_CS3(i)/N_users3;
    Mean_SINR_CS3(i) = Avg_SINR_CS3(i); 
    Mean_DD_CS3(i) = Avg_DD_CS3(i);
    Mean_CS_time3(i) = CS_time3(i);

    Err_Assoc_CS3(i) = std(Assoc_CS3(i))/N_users3;
    Err_SINR_CS3(i) = std(Avg_SINR_CS3(i)); 
    Err_DD_CS3(i) = std(Avg_DD_CS3(i));
    Err_CS_time3(i) = std(CS_time3(i));

end



%% Plots

XLabelsData = [{'Scenario 1'} ,{'Scenario 2'},{'Scenario 3'}];
colors = [0, 114, 178; 240, 228, 66; 213, 94, 0; 86, 180, 233; 0, 158, 115; 204, 121, 167]./255;

figure
hold on

colororder(colors)

ErrData1 = [Err_Assoc1(1), Err_Assoc_CS1(1), Err_Assoc1(2), Err_Assoc_CS1(2), Err_Assoc1(3), Err_Assoc_CS1(3);
    Err_Assoc2(1), Err_Assoc_CS2(1), Err_Assoc2(2), Err_Assoc_CS2(2), Err_Assoc2(3), Err_Assoc_CS2(3); ...
    Err_Assoc3(1), Err_Assoc_CS3(1), Err_Assoc3(2), Err_Assoc_CS3(2), Err_Assoc3(3), Err_Assoc_CS3(3)];

PlotData1 = [Mean_Assoc1(1), Mean_Assoc_CS1(1), Mean_Assoc1(2), Mean_Assoc_CS1(2), Mean_Assoc1(3), Mean_Assoc_CS1(3);
    Mean_Assoc2(1), Mean_Assoc_CS2(1), Mean_Assoc2(2), Mean_Assoc_CS2(2), Mean_Assoc2(3), Mean_Assoc_CS2(3);
    Mean_Assoc3(1), Mean_Assoc_CS3(1), Mean_Assoc3(2), Mean_Assoc_CS3(2), Mean_Assoc3(3), Mean_Assoc_CS3(3)];

bar_Assoc_Simulation = barwitherr(ErrData1,PlotData1);

set(gca,'XTick',[1, 2, 3])
set(gca,'XTickLabel',XLabelsData)
set(gca,"YTick", [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1])
set(gca,'YTickLabel',{'60%','65%','70%','75%','80%','85%','90%','95%','100%'})
legend('AHP (SINR) ','CS (SINR) ','AHP (TP) ','CS (TP) ','AHP (Coverage) ','CS (Coverage) ','Location', ...
    'southoutside','Orientation','Horizontal','Interpreter','latex')
ylim([0.6 1.02])
ylabel('Mean Coverage of Users')
yscale log
hold off

figure
hold on

colororder(colors)

ErrData2 = [Err_AHP_time1(1), Err_CS_time1(1), Err_AHP_time1(2), Err_CS_time1(2), Err_AHP_time1(3), Err_CS_time1(3);
    Err_AHP_time2(1), Err_CS_time2(1), Err_AHP_time2(2), Err_CS_time2(2), Err_AHP_time2(3), Err_CS_time2(3); ...
    Err_AHP_time3(1), Err_CS_time3(1), Err_AHP_time3(2), Err_CS_time3(2), Err_AHP_time3(3), Err_CS_time3(3)];

PlotData2 = [Mean_AHP_time1(1), Mean_CS_time1(1), Mean_AHP_time1(2), Mean_CS_time1(2), Mean_AHP_time1(3), Mean_CS_time1(3);
    Mean_AHP_time2(1), Mean_CS_time2(1), Mean_AHP_time2(2), Mean_CS_time2(2), Mean_AHP_time2(3), Mean_CS_time2(3);
    Mean_AHP_time3(1), Mean_CS_time3(1), Mean_AHP_time3(2), Mean_CS_time3(2), Mean_AHP_time3(3), Mean_CS_time3(3)];

% ErrData2 = [Err_AHP_time1(1), Err_AHP_time1(2), Err_AHP_time1(3),...
%     Err_CS_time1(1), Err_CS_time1(2), Err_CS_time1(3); ...
%     Err_AHP_time2(1), Err_AHP_time2(2), Err_AHP_time2(3), ...
%     Err_CS_time2(1), Err_CS_time2(2), Err_CS_time2(3); ...
%     Err_AHP_time3(1), Err_AHP_time3(2), Err_AHP_time3(3), ...
%     Err_CS_time3(1), Err_CS_time3(2), Err_CS_time3(3)];
% 
% PlotData2 = [Mean_AHP_time1(1), Mean_AHP_time1(2), Mean_AHP_time1(3), ...
%     Mean_CS_time1(1), Mean_CS_time1(2), Mean_CS_time1(3) ; ...
%     Mean_AHP_time2(1), Mean_AHP_time2(2), Mean_AHP_time2(3), ...
%     Mean_CS_time2(1), Mean_CS_time2(2), Mean_CS_time2(3);...
%     Mean_AHP_time3(1), Mean_AHP_time3(2), Mean_AHP_time3(3), ...
%     Mean_CS_time3(1), Mean_CS_time3(2), Mean_CS_time3(3)];

bar_Time_Simulation = barwitherr(ErrData2, PlotData2);

% ticks = round([0, max(Mean_AHP_time1)/4, ...
%     max(Mean_AHP_time1)/2, ...
%     3*max(Mean_AHP_time1)/4, ...
%     max(Mean_AHP_time1)], 1);

set(gca,'XTick',[1, 2, 3])
legend('AHP (SINR) ','CS (SINR) ','AHP (TP) ','CS (TP) ','AHP (Coverage) ','CS (Coverage) ','Location', ...
    'southoutside','Orientation','Horizontal','Interpreter','latex')
set(gca,'TickLength',[0 0]) 
set(gca,'XTickLabel',XLabelsData)
% set(gca,"YTick", ticks)
ylim([10^-2 10^5])
ylabel('Mean Running Time (sec)')
yscale log
hold off

figure
hold on

colororder(colors)

ErrData3 = [Err_SINR_AHP1(1), Err_SINR_CS1(1), Err_SINR_AHP1(2), Err_SINR_CS1(2), Err_SINR_AHP1(3), Err_SINR_CS1(3);
    Err_SINR_AHP2(1), Err_SINR_CS2(1), Err_SINR_AHP2(2), Err_SINR_CS2(2), Err_SINR_AHP2(3), Err_SINR_CS2(3); ...
    Err_SINR_AHP3(1), Err_SINR_CS3(1), Err_SINR_AHP3(2), Err_SINR_CS3(2), Err_SINR_AHP3(3), Err_SINR_CS3(3)];

PlotData3 = [Mean_SINR_AHP1(1), Mean_SINR_CS1(1), Mean_SINR_AHP1(2), Mean_SINR_CS1(2), Mean_SINR_AHP1(3), Mean_SINR_CS1(3);
    Mean_SINR_AHP2(1), Mean_SINR_CS2(1), Mean_SINR_AHP2(2), Mean_SINR_CS2(2), Mean_SINR_AHP2(3), Mean_SINR_CS2(3);
    Mean_SINR_AHP3(1), Mean_SINR_CS3(1), Mean_SINR_AHP3(2), Mean_SINR_CS3(2), Mean_SINR_AHP3(3), Mean_SINR_CS3(3)];

% ErrData3 = [Err_SINR_AHP1(1), Err_SINR_AHP1(2), Err_SINR_AHP1(3),...
%     Err_SINR_CS1(1), Err_SINR_CS1(2), Err_SINR_CS1(3); ...
%     Err_SINR_AHP2(1), Err_SINR_AHP2(2), Err_SINR_AHP2(3), ...
%     Err_SINR_CS2(1), Err_SINR_CS2(2), Err_SINR_CS2(3); ...
%     Err_SINR_AHP3(1), Err_SINR_AHP3(2), Err_SINR_AHP3(3), ...
%     Err_SINR_CS3(1), Err_SINR_CS3(2), Err_SINR_CS3(3)];
% 
% PlotData3 = [Mean_SINR_AHP1(1), Mean_SINR_AHP1(2), Mean_SINR_AHP1(3), ...
%     Mean_SINR_CS1(1), Mean_SINR_CS1(2), Mean_SINR_CS1(3); ...
%     Mean_SINR_AHP2(1), Mean_SINR_AHP2(2), Mean_SINR_AHP2(3), ...
%     Mean_SINR_CS2(1), Mean_SINR_CS2(2), Mean_SINR_CS2(3); ...
%     Mean_SINR_AHP3(1), Mean_SINR_AHP3(2), Mean_SINR_AHP3(3), ...
%     Mean_SINR_CS3(1), Mean_SINR_CS3(2), Mean_SINR_CS3(3)];

% PlotData3 = log2(1+(10.^((PlotData3)./10)));

bar_SINR_Simulation = barwitherr(ErrData3,PlotData3);

set(gca,'XTick',[1, 2, 3])
set(gca,'XTickLabel',XLabelsData)
% set(gca,"YTick", [0.9, 0.94, 0.96, 0.98, 1])
% set(gca,'YTickLabel',{'90%','94%','96%','98%','100%'})
legend('AHP (SINR) ','CS (SINR) ','AHP (TP) ','CS (TP) ','AHP (Coverage) ','CS (Coverage) ','Location', ...
    'southoutside','Orientation','Horizontal','Interpreter','latex')
% ylim([0 round(max(PlotData3,[],'all')+0.2)])
ylim([0 10])
ylabel('Mean SINR (dB)')
hold off

figure
hold on

colororder(colors)

ErrData4 = [Err_DD_AHP1(1), Err_DD_CS1(1), Err_DD_AHP1(2), Err_DD_CS1(2), Err_DD_AHP1(3), Err_DD_CS1(3);
    Err_DD_AHP2(1), Err_DD_CS2(1), Err_DD_AHP2(2), Err_DD_CS2(2), Err_DD_AHP2(3), Err_DD_CS2(3); ...
    Err_DD_AHP3(1), Err_DD_CS3(1), Err_DD_AHP3(2), Err_DD_CS3(2), Err_DD_AHP3(3), Err_DD_CS3(3)];

PlotData4 = [Mean_DD_AHP1(1), Mean_DD_CS1(1), Mean_DD_AHP1(2), Mean_DD_CS1(2), Mean_DD_AHP1(3), Mean_DD_CS1(3);
    Mean_DD_AHP2(1), Mean_DD_CS2(1), Mean_DD_AHP2(2), Mean_DD_CS2(2), Mean_DD_AHP2(3), Mean_DD_CS2(3);
    Mean_DD_AHP3(1), Mean_DD_CS3(1), Mean_DD_AHP3(2), Mean_DD_CS3(2), Mean_DD_AHP3(3), Mean_DD_CS3(3)];

% ErrData4 = [Err_DD_AHP1(1), Err_DD_AHP1(2), Err_DD_AHP1(3),...
%     Err_DD_CS1(1), Err_DD_CS1(2), Err_DD_CS1(3); ...
%     Err_DD_AHP2(1), Err_DD_AHP2(2), Err_DD_AHP2(3), ...
%     Err_DD_CS2(1), Err_DD_CS2(2), Err_DD_CS2(3); ...
%     Err_DD_AHP3(1), Err_DD_AHP3(2), Err_DD_AHP3(3), ...
%     Err_DD_CS3(1), Err_DD_CS3(2), Err_DD_CS3(3)];
% 
% PlotData4 = [Mean_DD_AHP1(1), Mean_DD_AHP1(2), Mean_DD_AHP1(3), ...
%     Mean_DD_CS1(1), Mean_DD_CS1(2), Mean_DD_CS1(3); ...
%     Mean_DD_AHP2(1), Mean_DD_AHP2(2), Mean_DD_AHP2(3), ...
%     Mean_DD_CS2(1), Mean_DD_CS2(2), Mean_DD_CS2(3); ...
%     Mean_DD_AHP3(1), Mean_DD_AHP3(2), Mean_DD_AHP3(3), ...
%     Mean_DD_CS3(1), Mean_DD_CS3(2), Mean_DD_CS3(3)];

% PlotData3 = log2(1+(10.^((PlotData3)./10)));

bar_DD_Simulation = barwitherr(ErrData4,PlotData4);

set(gca,'XTick',[1, 2, 3])
set(gca,'XTickLabel',XLabelsData)
% set(gca,"YTick", [0.9, 0.94, 0.96, 0.98, 1])
% set(gca,'YTickLabel',{'90%','94%','96%','98%','100%'})
legend('AHP (SINR) ','CS (SINR) ','AHP (TP) ','CS (TP) ','AHP (Coverage) ','CS (Coverage) ','Location', ...
    'southoutside','Orientation','Horizontal','Interpreter','latex')
% ylim([0 round(max(PlotData3,[],'all')+0.2)])
ylabel('Throughput, eMBB slice (Mbits/s)')
ylim([0 700])
hold off

Plot_time = toc;