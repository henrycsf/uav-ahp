clc
clear

s = rng;

N_users = 35000;
dist = 0.003;

X = 1270;

Y = 400;


vis = false;

Users = random_min_spacing(N_users, dist, vis) .* [X, Y];

Mu = randi(3,[N_users,1]);

for i = 1:1:length(Mu)
    
    if Mu(i) == 1
        
        Demand(i,1) = 6.0*10^6;
    end
    
    if Mu(i) == 2
        
        Demand(i,1) = 3.0*10^6;
    end

    if Mu(i) == 3
    
        Demand(i,1) = 0.8*10^6;
    end
    
end

dir = 'C:\Users\Noel Cafange\Documents\WORK\PeerJ_2\Trab final\';

save(horzcat(dir,'Users_festival.txt'), 'Users', '-ASCII');
save(horzcat(dir,'Demand_festival.txt'), 'Demand', '-ASCII');

figure(1)
hold on
scatter(Users(:,1), Users(:,2), '.')
rectangle('Position', [0 0 X Y])
axis([-50 X+50 -25 Y+25])
hold off