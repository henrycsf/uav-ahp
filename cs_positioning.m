function [bestnest, fmin, final_metrics, Simulation_time_CS]=cs_positioning(CaseSelect, Pt, N_users, X, Y, h, DataRate, ...
    Users, Demand, Zcov, Zcap, Rcov, Rcap, Nmi, Tol, n, alpha, beta, normal_values, plots)

tic

N_users = length(Users);

R_limit = min(Rcov, Rcap)*10^3;

UAV_limit = max(Zcov, Zcap);

UAV_Lb = [0 0];
UAV_Ub = [X Y];

N_UAV = UAV_limit;

Lb = [];
Ub = [];

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

nest = [];
Mean_TP_eMBB = [];
Mean_DD = [];
vetorbest = [];
vetormediafmin = [];
final_metrics = zeros(1,4);

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
while ((c < Nmi) && (- fmin < Tol || final_metrics(4) < normal_values(4))) || (c <= 2)
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
        bestnest=best;
        [fmin, SINR_max, final_metrics] = fobj(best);

    else
        if -fmin > Tol
            [fmin, SINR_max, final_metrics] = fobj(bestnest);
        end

    end
    
    vetorbest(c) = -fmin;
    vetormediafmin(c) = mean(-fitness);

    c = c + 1;
end %% End of iterations

Simulation_time_CS = toc;

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
    scatter(Users(:,1),Users(:,2), 40, SINR_max, "filled", "MarkerFaceAlpha", 0.85)
    xlabel('x (km)') 
    ylabel('y (km)')
    viscircles(CEN_CS,Rad, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'red');
    scatter(X_CS(1,:),Y_CS(1,:), 400, 'pentagram', 'MarkerEdgeColor', 'red', 'linewidth', 1.5)
    rectangle('Position', [0 0 X Y])
    axis([-X/10 X+X/10 -Y/10 Y+Y/10])
    daspect([1 1 1])
    xlabel('x (m)') 
    ylabel('y (m)')
    hold off

end

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
    [fnew, ~, ~] = fobj(newnest(j,:));
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
    function [z, SINR_max, metrics] = fobj(u)

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
TP_eMBB = 0;
eMBB_count = 0;

% for j = 1:1:size(SINR, 3)
% 
%     SINR_positioning(:,j) = SINR(:,j);
%     TP_positioning(:,j) = TP(:,j);
%     R_positioning(:,j) = R(:,j);
% 
% end

    Demand_max = max(Demand_Density,[],2);
    SINR_lin_max = max(SINR_lin,[],2);
    SINR_max = 10.*log10(SINR_lin_max);
    TP_max = max(TP,[],2);
    R_min = min(R,[],2);
    
    for k = 1:1:size(Users,1)
    
         if SINR_max(k,1) >= -10 && TP_max(k,1) > 0 && R_min(k,1) <= R_limit
            
             Assoc = Assoc + 1;
            
            if User_Service(k) == 1
                TP_eMBB = TP_eMBB + TP_max(k);
                 eMBB_count = eMBB_count + 1;
            end

        end
    
    end
    
    Mean_SINR_lin = mean(SINR_lin_max,1);
    
    Mean_SINR = 10*log10(Mean_SINR_lin);

    Mean_TP_eMBB = TP_eMBB/eMBB_count;
    
    Mean_DD = mean(Demand_max,1);
    
    % Cost Functions
    
    z1 = -Mean_SINR/normal_values(1);
    
    z2 = -Mean_TP_eMBB/normal_values(2);
    
    z3 = -Mean_DD/normal_values(3);
    
    z4 = -Assoc/normal_values(4);

    switch CaseSelect
        case 1
            z = 0.5559*z1 + 0.1364*z2 + 0.0489*z3 + 0.2589*z4;
        case 2
            z = 0.0834*z1 + 0.6259*z2 + 0.2229*z3 + 0.0678*z4;
        case 3
            z = 0.1558*z1 + 0.0856*z2 + 0.0491*z3 + 0.7095*z4;
        otherwise
            error("CaseSelect value not allowed. Please put 1, 2 or 3.")
    end

    metrics = [Mean_SINR, Mean_TP_eMBB, Mean_DD, Assoc];
end

end