clear
clc

nmi = 1000;

Tol = 1;

plots = false;

% CaseSelect:

% 1 for SINR prio
% 2 for TP prio
% 3 for users with higher demand prio

% This will give objective functions based on the criteria weights
% given by AHP in each of those 3 options.

CaseSelect = 1;

dir = '/home/henry-ferreira/Documentos/WORK/PeerJ_2/Trab final/';

% % Adaptive-AHP and CS, Simulation 2

DataRate2 = [100*10^6, 15*10^6, 7*10^6];

Users2 = {load(horzcat(dir,'Users2_02.txt')), load(horzcat(dir,'Users2_01.txt')), load(horzcat(dir,'Users2_0.txt')), load(horzcat(dir,'Users2.txt')), ...
    load(horzcat(dir,'Users2_2.txt')), load(horzcat(dir,'Users2_3.txt')), load(horzcat(dir,'Users2_4.txt')), load(horzcat(dir,'Users2_5.txt'))};

Demand2 = {load(horzcat(dir,'Demand2_02.txt')),load(horzcat(dir,'Demand2_01.txt')),load(horzcat(dir,'Demand2_0.txt')), load(horzcat(dir,'Demand2.txt')), ...
    load(horzcat(dir,'Demand2_2.txt')), load(horzcat(dir,'Demand2_3.txt')), load(horzcat(dir,'Demand2_4.txt')), load(horzcat(dir,'Demand2_5.txt'))};

delta_cov2 = 0.40;
X2 = 200;
Y2 = 200;
h2 = 30; %height in meters

alpha2 = 27.23;
beta2 = 0.08;

for k = 1:length(Users2)

N_users2(k) = length(Users2{1,k});

[Zcov2(k),Zcap2(k),Rcov2(k),Rcap2(k),Pt2(k)] = numerology(h2,X2,Y2,N_users2(k),delta_cov2, DataRate2, alpha2, beta2);

[UAV2{1,k}, Best2{1,k}, AHP_time2{1,k}, Assoc2{1,k}, Mean_SINR_AHP2{1,k}, Mean_DD_AHP2{1,k}] = drone_positioning_AHP(CaseSelect, Pt2(k), N_users2(k), X2, ...
    Y2, h2, DataRate2, Users2{1,k}, Demand2{1,k}, Zcov2(k), Zcap2(k), Rcov2(k), Rcap2(k), nmi, alpha2, beta2, plots);

[bestnest2{1,k}, fmin2{1,k}, Assoc_CS2{1,k}, Mean_SINR_CS2{1,k}, Mean_DD_CS2{1,k}, Simulation_time_CS2{1,k}]=cs_positioning(CaseSelect, Pt2(k), N_users2(k), X2, ...
    Y2, h2, DataRate2, Users2{1,k}, Demand2{1,k}, Zcov2(k), Zcap2(k), Rcov2(k), Rcap2(k), nmi, Tol, alpha2, beta2, plots);

AHP_time(k) = AHP_time2{1,k};
CS_time(k) = Simulation_time_CS2{1,k};

Assoc_AHP(k) = Assoc2{1,k};
Assoc_CS(k) = Assoc_CS2{1,k};

SINR_AHP(k) = Mean_SINR_AHP2{1,k};
SINR_CS(k) = Mean_SINR_CS2{1,k};

end

figure
hold on

b = bar(N_users2,[Zcap2; Pt2],'stacked','LineWidth',0.1,'FaceAlpha',0.6);
b(1).FaceColor = 'blue';
b(2).FaceColor = 'cyan';

yyaxis left
plot(N_users2, 10^3.*Rcap2,'k','LineWidth',2.0)
%xlabel("Number of Users (UEs)")
xticks([500, 1000, 1500, 2000, 2500, 3000, 3500, 4000])
xticklabels(["Z_{lim} = 3 " + "\newline" + "P_{t} = 23 dB" + "\newline" + "N_{UE} = 500", ...
    "Z_{lim} = 5" + "\newline" + "P_{t} = 20dB" + "\newline" + "N_{UE} = 1000", ...
    "Z_{lim} = 8" + "\newline" + "P_{t} = 18dB" + "\newline" + "N_{UE} = 1500", ...
    "Z_{lim} = 10" + "\newline" + "P_{t} = 17dB" + "\newline" + "N_{UE} = 2000", ...
    "Z_{lim} = 13" + "\newline" + "P_{t} = 17dB" + "\newline" + "N_{UE} = 2500", ...
    "Z_{lim} = 15" + "\newline" + "P_{t} = 16dB" + "\newline" + "N_{UE} = 3000", ...
    "Z_{lim} = 18" + "\newline" + "P_{t} = 15dB" + "\newline" + "N_{UE} = 3500", ...
    "Z_{lim} = 20" + "\newline" + "P_{t} = 15dB" + "\newline" + "N_{UE} = 4000"]);
ylabel('Cell Limit Radius (m)')

yyaxis right
plot(N_users2, AHP_time,'r','LineWidth',2.0)
plot(N_users2, CS_time, 'r','LineWidth',2.0,'LineStyle','--')
ylabel('Simulation Time (sec.)','Color','black')
yscale("log")

legend(["$Z_{lim}$","$P_{t}$","Cell Radius $R_{lim}$","Running Time (UAV-AHP)","Running Time (CS)"],'FontSize',12,'Interpreter','latex','Location', ...
    'southoutside','Orientation','Horizontal')

a.TickLabelInterpreter = "latex";
xlim([250, 4250])

hold off

figure
hold on 

yyaxis left
plot(N_users2, SINR_AHP,'LineWidth',1.5)
plot(N_users2, SINR_CS,'LineWidth',1.5,'LineStyle','--')
xlabel("Number of Users (UEs)")
ylabel('Average SINR (dB)')

yyaxis right
plot(N_users2, 1-(Assoc_AHP./N_users2),'LineWidth',1.5)
plot(N_users2, 1-(Assoc_CS./N_users2),'LineWidth',1.5,'LineStyle','--')
ylabel('User Outage')

legend(["SINR (UAV-AHP)","SINR (CS)","Outage (UAV-AHP)","Outage (CS)"],'FontSize',12,'Interpreter','latex')

a.TickLabelInterpreter = "latex";

hold off