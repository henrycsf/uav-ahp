clc
clear

s = rng;

N_users = 500;
dist = 0.02;

XY = 200;

UserXY = XY - XY/10;

vis = false;

Users = XY/20 + UserXY * random_min_spacing(N_users, dist, vis);

Mu = randi(3,[N_users,1]);

for i = 1:1:length(Mu)
    
    if Mu(i) == 1
        
        Demand(i,1) = 100*10^6;
    end
    
    if Mu(i) == 2
        
        Demand(i,1) = 15*10^6;
    end
    
    if Mu(i) == 3
        
        Demand(i,1) = 7*10^6;
    end
    
end

dir = '/home/henry-ferreira/Documentos/WORK/PeerJ_2/Trab final/';

save(horzcat(dir,'Users2_02.txt'), 'Users', '-ASCII');
save(horzcat(dir,'Demand2_02.txt'), 'Demand', '-ASCII');

figure(1)
hold on
rectangle('Position', [0 0 XY XY])
scatter(Users(:,1), Users(:,2))
axis([0 XY 0 XY])
hold off