clear
clc

dir = '/home/henry-ferreira/Documentos/WORK/PeerJ_2/Trab final/';

CaseSelect = 1;

nmi = 100;

delta = 0.05 : 0.05 : 0.95;

DataRate2 = [100*10^6, 15*10^6, 7*10^6];

Users2 = load(horzcat(dir,'Users2.txt'));

Demand2 = load(horzcat(dir,'Demand2.txt'));

plots = false;

N_users2 = length(Users2);
X2 = 200;
Y2 = 200;
h2 = 40; %height in meters

DataPoints = 10;

for i = 1:length(delta)
    for j = 1:DataPoints
    
    [Zcov2(i,j),Zcap2(i,j),Rcov2(i,j),Rcap2(i,j), Pt2(i,j)] = numerology(h2,X2,Y2,N_users2,delta(i), DataRate2);
    
    [~, ~, AHP_time2(i,j), Assoc2(i,j), Mean_SINR_AHP2(i,j), Mean_SINR_1(i,j), ...
        Mean_SINR_2(i,j), Mean_SINR_3(i,j)] = drone_positioning_AHP_loads(CaseSelect, Pt2(i,j), N_users2, X2, ...
        Y2, h2, DataRate2, Users2, Demand2, Zcov2(i,j), Zcap2(i,j), Rcov2(i,j), Rcap2(i,j), nmi, delta(i), plots);
    
    end
end

Mean_SINR = mean(Mean_SINR_AHP2,2);
Mean_SINR_eMBB = mean(Mean_SINR_1,2);
Mean_SINR_mMTC = mean(Mean_SINR_2,2);
Mean_SINR_URLLC = mean(Mean_SINR_3,2);

Mean_SE = log2(1+(10.^((Mean_SINR)./10)));
Mean_SE_eMBB = log2(1+(10.^((Mean_SINR_eMBB)./10)));
Mean_SE_mMTC = log2(1+(10.^((Mean_SINR_mMTC)./10)));
Mean_SE_URLLC = log2(1+(10.^((Mean_SINR_URLLC)./10)));

figure
hold on
plot(delta, Mean_SE)
plot(delta, Mean_SE_eMBB)
plot(delta, Mean_SE_mMTC)
plot(delta, Mean_SE_URLLC)
xlim([delta(1) delta(i)])
ylim([0 max(Mean_SE)+0.05])
hold off