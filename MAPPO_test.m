clc;
clear;
close all;

%% Load Trained Agents and Critic
load("MAPPO_UAV_Agents.mat", "agents", "critic");

%% Environment Setup
load("numerology_festival.mat");

numUAVs = max(Zcov, Zcap); 
Rlimit = min(Rcov, Rcap)*10^3;
stepSize = 0.3 * [ areaSizeX/max([areaSizeX areaSizeY])*Rlimit, ...
                     areaSizeY/max([areaSizeX areaSizeY])*Rlimit ];

% Create the environment using the new MAPPO_Environment class.
env = MAPPO_Environment(numUAVs, numUEs, UEPositions, demandVector, height, areaSizeX, areaSizeY, stepSize);

%% Validation Parameters
numTestEpisodes = 20;
maxTimeSteps = env.maxSteps;  % e.g., 9 time steps per episode

episodeRewards = zeros(numTestEpisodes,1);
coverageHistory = zeros(numTestEpisodes,1);
SINRHistory = zeros(numTestEpisodes,1);
TPHistory = zeros(numTestEpisodes,1);
demandHistory = zeros(numTestEpisodes,1);

for ep = 1:numTestEpisodes
    % Reset the environment (global state vector [numUAVs*2,1])
    state = reset(env);
    
    % Log metrics at the start of the episode:
    [coverageStart, SINRStart, throughputStart, demandStart] = env.evaluateCoverage();
    fprintf('Episode %d Start: Coverage = %.2f, SINR = %.2f, Throughput = %.2f, Demand = %.2f\n',...
        ep, coverageStart, SINRStart, throughputStart, demandStart);
    
    epReward = 0;
    
    for t = 1:maxTimeSteps
        % For each UAV, obtain its action using its actor and aggregate them.
        globalAction = zeros(numUAVs*2,1);
        for i = 1:numUAVs
            % Each agent's actor outputs a global action vector.
            action_i = getAction(agents{i}, state);  % Expected to return a numeric array.
            globalAction = globalAction + action_i{1};
        end
        
        % Step the environment with the aggregated global action.
        [nextState, reward, done] = step(env, globalAction);
        epReward = epReward + reward;
        
        % Visualization: display UAV and UE positions.
        UAVPositions = reshape(nextState, numUAVs, 2);
        figure(1);
        clf;
        hold on;
        scatter(UAVPositions(:,1), UAVPositions(:,2), 100, 'bo', 'filled');
        scatter(env.UEPositions(:,1), env.UEPositions(:,2), 25, 'rx');
        xlim([0, env.areaSizeX]);
        ylim([0, env.areaSizeY]);
        title(sprintf('Episode %d, Time Step %d', ep, t));
        xlabel('X Position'); ylabel('Y Position');
        legend('UAVs', 'UEs');
        drawnow;
        pause(0.5);
        
        state = nextState;
        if done
            fprintf('Episode %d finished at time step %d.\n', ep, t);
            break;
        end
    end
    
    % Log metrics at the end of the episode:
    [coverageEnd, SINREnd, throughputEnd, demandEnd] = env.evaluateCoverage();
    fprintf('Episode %d End: Coverage = %.2f, SINR = %.2f, Throughput = %.2f, Demand = %.2f\n',...
        ep, coverageEnd, SINREnd, throughputEnd, demandEnd);
    
    episodeRewards(ep) = epReward;
    coverageHistory(ep) = coverageEnd;
    SINRHistory(ep) = SINREnd;
    TPHistory(ep) = throughputEnd;
    demandHistory(ep) = demandEnd;
    
    fprintf('Episode %d: Total Reward = %.2f\n', ep, epReward);
end

%% Plot Overall Episode Metrics
figure(2);
subplot(3,1,1);
plot(1:numTestEpisodes, episodeRewards, '-o','LineWidth',1.5);
xlabel('Episode');
ylabel('Total Reward');
title('Episode Rewards');

subplot(3,1,2);
plot(1:numTestEpisodes, coverageHistory, '-o','LineWidth',1.5);
xlabel('Episode');
ylabel('Coverage');
title('Episode Coverage');

subplot(3,1,3);
plot(1:numTestEpisodes, SINRHistory, '-o','LineWidth',1.5);
xlabel('Episode');
ylabel('SINR');
title('Episode SINR');

figure(3);
subplot(1,2,1);
plot(1:numTestEpisodes, TPHistory, '-o','LineWidth',1.5);
xlabel('Episode');
ylabel('Throughput');
title('Episode Throughput');

subplot(1,2,2);
plot(1:numTestEpisodes, demandHistory, '-o','LineWidth',1.5);
xlabel('Episode');
ylabel('Demand');
title('Episode Demand');
