function [bestnest,fmin, Assoc, Mean_SINR, Mean_TP_eMBB, Simulation_time_CS]=cs_positioning(CaseSelect, Pt, N_users, X, Y, h, DataRate, ...
    Users, Demand, Zcov, Zcap, Rcov, Rcap, Nmi, Tol, alpha, beta, plots)

tic

N_users = length(Users);

R_limit = min(Rcov, Rcap)*10^3;

UAV_limit = max(Zcov, Zcap);

UAV_Lb = [0 0];
UAV_Ub = [X Y];

N_UAV = UAV_limit;

Lb = [];
Ub = [];

n = 25;

% Discovery rate of alien eggs/solutions
pa = 0.25;
%D;T1;T2;W1;W2

%% Simple bounds of the search domain

for r = 1:1:N_UAV
% Lower bounds
Lb = horzcat(Lb,UAV_Lb);
% Upper bounds
Ub = horzcat(Ub,UAV_Ub);
end

nest = ones(n, length(Lb));

% Random initial solutions
for i=1:n
    nest(i, :)=round(Lb+(Ub-Lb).*rand(size(Lb)));
end

% Get the current best
fitness=10^10*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness);
N_iter=0;

%% Starting iterations
c = 1;
while (- fmin < Tol && c < Nmi)
    % Generate new solutions (but keep the current best)
     new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
     [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
     
    % Update the counter
      N_iter=N_iter+n; 
      
    % Discovery and randomization
      new_nest=empty_nests(nest,Lb,Ub,pa);
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
      
    % Update the counter again
      N_iter=N_iter+n;
      
    % Find the best objective so far  
    if fnew<fmin
        fmin=fnew;
        bestnest=best;

        [fnew, Assoc, Mean_SINR, SINR_max] = fobj(best);

    end
    
    vetorbest(c) = -fmin;
    vetormediafmin(c) = mean(-fitness);

    c = c + 1;
end %% End of iterations

%% Post-optimization processing
%% Display all the nests

%%%%%%%% PLOT %%%%%%%%%%

if plots == true

    figure
    plot(1:c-1,vetorbest,'k',1:c-1,vetormediafmin,'b--');
    legend({'Best Fitness Value','Average Fitness Value'}, 'FontSize', 12)
    xlabel('Iterations')
    ylabel('Fitness')
    
    X_CS = [];
    Y_CS = [];
    
    for i = 1:2:length(bestnest)
    
        X_CS = horzcat(X_CS,bestnest(i));
        Y_CS = horzcat(Y_CS,bestnest(i+1));
    
    end
    
    CEN_CS = transpose([X_CS; Y_CS]);
    
    for i = 1:1:length(bestnest)/2
        Rad(1,i) = R_limit;
    end
    
    figure
    hold on
    colorbar
    title('SINR (dB) by Cuckoo Search');
    scatter(Users(:,1),Users(:,2), 20, SINR_max, "filled")
    xlabel('x (km)') 
    ylabel('y (km)')
    viscircles(CEN_CS,Rad, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'red');
    scatter(X_CS(1,:),Y_CS(1,:), 250, 'diamond', 'MarkerEdgeColor', 'red', 'linewidth', 2)
    rectangle('Position', [0 0 X Y])
    axis([-X/10 X+X/10 -Y/10 Y+Y/10])
    daspect([1 1 1])
    xlabel('x (m)') 
    ylabel('y (m)')
    hold off

end

Simulation_time_CS = toc;

%% --------------- All subfunctions are listed below ------------------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);

% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).

Beta=1/2;
sigma=(gamma(1+Beta)*sin(pi*Beta/2)/(gamma((1+Beta)/2)*Beta*2^((Beta-1)/2)))^(1/Beta);

for j=1:n
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=round(randn(size(s))*sigma);
    v=randn(size(s));
    step=u./abs(v).^(1/Beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale; 
    % otherwise, Levy flights may become too aggresive/efficient, 
    % which makes new solutions (even) jump out side of the design domain 
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
   % Apply simple bounds/limits
   nest(j,:)=simplebounds(s,Lb,Ub);
end

end

%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness)
% Evaluating all new solutions
for j=1:size(nest,1)
    [fnew, Assoc, Mean_SINR, SINR_max] = fobj(newnest(j,:));
    if fnew<=fitness(j)
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness);
best=nest(K,:);

end

%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;

% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions.  Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes.  
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=round(nest+stepsize.*K);
for j=1:size(new_nest,1)
    s=new_nest(j,:);
  new_nest(j,:)=simplebounds(s,Lb,Ub);  
end

end

% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
  
end

%% You can replace the following by your own functions
% A d-dimensional objective function
function [z, Assoc, Mean_SINR, SINR_max] = fobj(u)

DronePop = zeros(2,N_UAV);

for e = 1:2:2*N_UAV
    DronePop(:,(e+1)/2) = [u(1,e);u(1,e+1)];
end

N_users = length(Users);

f = 3.5 * 10^9;
velc = 299792458;
ZetaLOS = 1;
ZetaNLOS = 20;
B = [50*10^6 20*10^6 10*10^6];
q = -174 + 10*log10(B);

Gt = 3;
Gr = 0;

for b = 1:1:N_UAV
    for k = 1:1:N_users

        D(k,b) = sqrt(((Users(k,1)-DronePop(1,b)).^2) + (Users(k,2)-DronePop(2,(b))).^2 + h.^2);
        R(k,b) = sqrt(abs(((Users(k,1)-DronePop(1,b)).^2) + (Users(k,2)-DronePop(2,(b))).^2));
        Demand_Density(k,b) = Demand(k)./(D(k,b).^2);

         switch Demand(k)
    
             case DataRate(1)
                 User_Service(k) = 1;
                 
             case DataRate(2)
                 User_Service(k) = 2;
    
             case DataRate(3)
                 User_Service(k) = 3;
    
             otherwise
                 error("Out of service");
         end

    end
end

    theta = atan(h./R);

    Z = (alpha*exp(-beta*((180/pi).*theta - alpha)));

    PL = 20*log10((4*pi*f*D./velc))+((ZetaLOS+Z.*ZetaNLOS)./(1+Z));

    Pr_user = Pt - PL + Gt + Gr;

    Pr_Linear = (10.^((Pr_user-30)./10));

    % SINR Estimation

    for k = 1:1:N_UAV

        for j = 1:1:N_users

            inter = 0;

                for m = 1:1:N_UAV

                    if (k ~= m)
    
                        inter = inter + Pr_Linear(j,m);
    
                    end
                end

            Interference_lin(j,k) = inter;

            Interference(j,k) = 10*log10(inter);

            switch Demand(j)

                case DataRate(1)
                if R(j,k) > R_limit
                    SINR_lin(j,k) = 0;
                    TP(j,k) = 0;
                else
                    SINR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(1)-30)/10)) + inter);
                    SNR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(1)-30)/10)));
                    TP(j,k) = (10^-6) * (B(1))*log2(1+SNR_lin(j,k));
                end

                SINR(j,k) = 10.*log10(SINR_lin(j,k));

                case DataRate(2)

                if R(j,k) > R_limit
                    SINR_lin(j,k) = 0;
                    TP(j,k) = 0;
                else
                    SINR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(2)-30)/10)) + inter);
                    SNR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(2)-30)/10)));
                    TP(j,k) = (10^-6) * (B(2))*log2(1+SNR_lin(j,k));
                end

                SINR(j,k) = 10.*log10(SINR_lin(j,k));

                case DataRate(3)

                if R(j,k) > R_limit
                    SINR_lin(j,k) = 0;
                    TP(j,k) = 0;
                else
                    SINR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(3)-30)/10)) + inter);
                    SNR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(3)-30)/10)));
                    TP(j,k) = (10^-6) * (B(3))*log2(1+SNR_lin(j,k));
                end

                SINR(j,k) = 10.*log10(SINR_lin(j,k));

            end

        end

    end

Assoc = 0;

% for j = 1:1:size(SINR, 3)
% 
%     SINR_positioning(:,j) = SINR(:,j);
%     TP_positioning(:,j) = TP(:,j);
%     R_positioning(:,j) = R(:,j);
% 
% end

Demand_max = max(Demand_Density,[],2);
SINR_lin_max = max(SINR_lin,[],2);
SINR_max = max(SINR,[],2);
TP_max = max(TP,[],2);
R_min = min(R,[],2);

for k = 1:1:size(Users,1)

     if SINR_max(k,1) >= -10 && TP_max(k,1) > 0 && R_min(k,1) <= R_limit
        
         Assoc = Assoc + 1;
    end

end

SINR_lin_normal = (SINR_lin - min(SINR_lin,[],2))./(max(SINR_lin,[],2)-min(SINR_lin,[],2));

TP_normal = (TP - min(TP,[],2))./(max(TP,[],2)-min(TP,[],2));

Demand_normal = (Demand_Density - min(Demand_Density,[],2))./(max(Demand_Density,[],2)-min(Demand_Density,[],2));

Coverage_Normal = (Assoc - 0) ./ (N_users - 0);

Mean_SINR_lin = mean(SINR_lin_max,1);

Mean_SINR = 10*log10(Mean_SINR_lin);

Mean_TP = mean(TP_max,1);

TP_eMBB = 0;
Assoc_eMBB = 0;

for k = 1:length(TP_max)
    if User_Service(k) == 1

        TP_eMBB = TP_eMBB + TP_max(k);
        Assoc_eMBB = Assoc_eMBB + 1;
    end
end

Mean_TP_eMBB = TP_eMBB/Assoc_eMBB;

Mean_DD = mean(Demand_max,1);

% Cost Functions

z1 = -mean(max(SINR_lin_normal,[],2),1,"omitnan");

z2 = -mean(max(TP_normal,[],2),1,"omitnan");

z3 = -mean(max(Demand_normal,[],2),1,"omitnan");

z4 = -mean(max(Coverage_Normal,[],2),1,"omitnan");

    switch CaseSelect
        case 1
            z = 0.5559*z1 + 0.1364*z2 + 0.0489*z3 + 0.2589*z4;
        case 2
            z = 0.1404*z1 + 0.5805*z2 + 0.1365*z3 + 0.1427*z4;
        case 3
            z = 0.1558*z1 + 0.0856*z2 + 0.0491*z3 + 0.7095*z4;
        otherwise
            error("CaseSelect value not allowed. Please put 1, 2 or 3.")
    end

end

end