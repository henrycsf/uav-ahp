clc;
clear;
close all;

%% Load the Trained Agent
load('MADQN_Agent.mat', 'agent');

%% Environment Setup (Same as during training)
load("numerology_festival.mat");

numUAVs = max(Zcov, Zcap); 
Rlimit = min(Rcov, Rcap)*10^3;

stepSize = 0.3*[areaSizeX/max([areaSizeX areaSizeY]) * Rlimit, areaSizeY/(max([areaSizeX areaSizeY])) * Rlimit];

[UAVPositions, ~, ~, ~] = drone_positioning_AHP(CaseSelect, Pt, numUEs, areaSizeX, ...
    areaSizeY, height, DataRate, UEPositions, demandVector, Zcov, Zcap, Rcov, Rcap, nmi, alpha, beta, plots);

% Create the environment (ensure MADQN_Environment.m is modified for 5 actions per UAV)
% For 4 UAVs with 5 actions each, the joint action space has 5^4 = 625 possible actions.
env = MADQN_Environment(numUAVs, numUEs, UEPositions, UAVPositions, demandVector, height, areaSizeX, areaSizeY, stepSize);

%% Validation Parameters
numTestEpisodes = 30;
maxTimeSteps = 9;  % renamed from "step" to avoid conflict

%obj.UEPositions = obj.UEPositions + (4*randn(obj.numUEs,2));

% Arrays to collect episode rewards, coverage, SINR, throughput, and demand
episodeRewards = zeros(numTestEpisodes,1);
coverageHistory = zeros(numTestEpisodes,1);
SINRHistory = zeros(numTestEpisodes,1);
TPHistory = zeros(numTestEpisodes,1);
demandHistory = zeros(numTestEpisodes,1);

for ep = 1:numTestEpisodes
    % Reset environment at the beginning of each episode
    state = reset(env);
    epReward = 0;
    
    fprintf('\n=== Episode %d ===\n', ep);
    
    for t = 1:maxTimeSteps
        % Select an action using the trained agent.
        % (For discrete DQN, getAction should return a discrete value.)
        actionCell = getAction(agent, state);
        action = actionCell{1};
        fprintf('Time Step %d: Action = %d\n', t, action);
        
        % Step the environment using the chosen action.
        [nextState, reward, done] = step(env, action);
        epReward = epReward + reward;
        
        % Extract UAV positions from the state vector.
        uavPositions = reshape(nextState, numUAVs, 2);
        
        % Visualization:
        figure(1);
        clf;
        hold on;
        % Plot UAV positions as blue circles.
        scatter(uavPositions(:,1), uavPositions(:,2), 100, 'bo', 'filled');
        % Plot UE positions as red crosses.
        scatter(env.UEPositions(:,1), env.UEPositions(:,2), 25, 'rx');
        xlim([0 areaSize]);
        ylim([0 areaSize]);
        title(sprintf('Episode %d, Time Step %d', ep, t));
        xlabel('X Position (m)'); ylabel('Y Position (m)');
        legend('UAVs', 'UEs');
        drawnow;
        hold off;
        
        % Update state for next time step.
        state = nextState;
        
        % If episode is finished, break.
        if done
            fprintf('Episode %d finished at time step %d.\n', ep, t);
            break;
        end
        
        pause(0.25); % Brief pause for visualization clarity.
    end
    
    % Evaluate performance metrics after each episode.
    [coverage, SINR, throughput, demand] = env.evaluateCoverage();
    episodeRewards(ep) = epReward;
    coverageHistory(ep) = coverage;
    SINRHistory(ep) = SINR;
    TPHistory(ep) = throughput;
    demandHistory(ep) = demand;
    
    fprintf('Episode %d: Total Reward = %.2f, Coverage = %.2f, SINR = %.2f, TP = %.2f, Demand = %.2f\n', ...
        ep, epReward, coverage, SINR, throughput, demand);
end

%% Plot Episode Statistics
figure(2);
subplot(3,1,1);
plot(1:numTestEpisodes, episodeRewards, '-o');
xlabel('Episode'); ylabel('Total Reward');
title('Episode Rewards');

subplot(3,1,2);
plot(1:numTestEpisodes, coverageHistory, '-o');
xlabel('Episode'); ylabel('Coverage');
title('Episode Coverage');

subplot(3,1,3);
plot(1:numTestEpisodes, SINRHistory, '-o');
xlabel('Episode'); ylabel('SINR');
title('Episode SINR');

figure(3);
subplot(1,2,1);
plot(1:numTestEpisodes, TPHistory, '-o');
xlabel('Episode'); ylabel('Throughput (TP)');
title('Episode Throughput');

subplot(1,2,2);
plot(1:numTestEpisodes, demandHistory, '-o');
xlabel('Episode'); ylabel('Demand Density');
title('Episode Demand Density');
