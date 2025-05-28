clc;
clear;
close all;

% pool = gcp('nocreate'); % Check if a parallel pool is already running
% if isempty(pool)
%     poolobj.NumWorkers = 12;
%     pool = parpool("Processes",poolobj.NumWorkers); % Start a pool with 8 workers
% end

%% Environment Setup
load("numerology_festival.mat");

numUAVs = max(Zcov, Zcap); 
cellRadius = min(Rcov, Rcap)*10^3;

stepSize = 0.3 * [areaSizeX/max([areaSizeX areaSizeY])*cellRadius, areaSizeY/max([areaSizeX areaSizeY])*cellRadius];

% Create the environment using the new MAPPO_Environment class.
env = MAPPO_Environment(numUAVs, numUEs, UEPositions, demandVector, height, areaSizeX, areaSizeY, stepSize);

%% Global Critic Setup (Centralized Critic)
% Global observation: concatenated UAV positions [numUAVs*2, 1]
globalObsDim = numUAVs * 2;  
upperLimit = repmat([areaSizeX; areaSizeY], numUAVs, 1);
globalObsInfo = rlNumericSpec([globalObsDim 1], 'LowerLimit', zeros(globalObsDim,1), 'UpperLimit', upperLimit);
globalObsInfo.Name = 'globalState';

% Critic Network (global state input)
criticLayers = [
    featureInputLayer(globalObsDim, "Normalization", "none", "Name", "globalState")
    fullyConnectedLayer(128, "Name", "crit_fc1")
    reluLayer("Name", "crit_relu1")
    fullyConnectedLayer(64, "Name", "crit_fc2")
    reluLayer("Name", "crit_relu2")
    fullyConnectedLayer(1, "Name", "value_out")
    ];

criticNet = layerGraph(criticLayers);
critic = rlValueFunction(criticNet, globalObsInfo, 'ObservationInputNames', {'globalState'});
critic.UseDevice = "gpu";

%% Actor Setup for Each UAV (Decentralized Actor with Local Extraction & Embedding)
% We want the actor to receive the global state (for compatibility) but then extract only its local [x;y]
% and then output a local action [dx;dy] that is embedded into a global action vector.
% The environment's action spec is global: [numUAVs*2,1].
%
% We'll define the actor network for each UAV using a helper function.
%
% The global observation spec is used for both actor and critic to satisfy MATLAB's requirements.
actorObsInfo = globalObsInfo;  
% Actor action info remains as a global action vector:
actorActInfo = rlNumericSpec([numUAVs*2, 1], 'LowerLimit', -repmat(stepSize(:), numUAVs, 1), 'UpperLimit', repmat(stepSize(:), numUAVs, 1));
actorActInfo.Name = 'actorAction';

%% Create a function to build an actor network for UAV i.
% This network starts with a global input layer (size [numUAVs*2,1]),
% then uses a LocalStateExtractionLayer (named "localExtract") to pick out UAV i's [x;y],
% then processes it through shared layers, then outputs a 2-D local action.
% Finally, a LocalActionEmbeddingLayer (named "global_mean" and "global_std") maps the local output
% into a global action vector (of dimension [numUAVs*2,1]) by inserting the local action in the correct indices.
createActorNetForAgent = @(agentIndex, numUAVs, stepSize) localActorNet(agentIndex, numUAVs, stepSize);

%% Local function to build an actor network for UAV agent
function net = localActorNet(agentIndex, numUAVs, stepSize)
    % Global input layer: matches global observation [numUAVs*2, 1]
    globalInput = featureInputLayer(numUAVs*2, "Normalization", "none", "Name", "globalState");
    
    % Local extraction: extracts the [2,1] portion for UAV 'agentIndex'
    localExtract = LocalStateExtractionLayer(agentIndex, "localExtract");
    
    % Shared layers after extraction (local state now of size [2,1])
    sharedLayers = [
        fullyConnectedLayer(64, "Name", "actor_fc1")
        reluLayer("Name", "actor_relu1")
        fullyConnectedLayer(32, "Name", "actor_fc2")
        reluLayer("Name", "actor_relu2")
        ];
    
    % Mean output branch (local action)
    meanBranch = [
        fullyConnectedLayer(2, "Name", "action_mean_fc")
        tanhLayer("Name", "tanh")
        scalingLayer("Name", "local_mean", "Scale", stepSize(:))];
    
    % Standard deviation branch (local action)
    stdBranch = [
        fullyConnectedLayer(2, "Name", "std_dev_fc")
        softplusLayer("Name", "local_std")];
    
    % Now embed the local outputs into a global action vector.
    % The embedding layer will take the local 2-D vector and produce a global vector of zeros except at positions corresponding to this UAV.
    embeddingMean = LocalActionEmbeddingLayer(agentIndex, numUAVs, "global_mean");
    embeddingStd = LocalActionEmbeddingLayer(agentIndex, numUAVs, "global_std");
    
    % Assemble the layer graph.
    net = layerGraph();
    net = addLayers(net, globalInput);
    net = addLayers(net, localExtract);
    net = addLayers(net, sharedLayers);
    net = addLayers(net, meanBranch);
    net = addLayers(net, stdBranch);
    net = addLayers(net, embeddingMean);
    net = addLayers(net, embeddingStd);
    
    % Connect global input to local extraction.
    net = connectLayers(net, "globalState", "localExtract");
    % Connect local extraction to shared layers.
    net = connectLayers(net, "localExtract", "actor_fc1");
    % Connect shared layers to both output branches.
    net = connectLayers(net, "actor_relu2", "action_mean_fc");
    net = connectLayers(net, "actor_relu2", "std_dev_fc");
    % Connect mean branch to embedding layer.
    net = connectLayers(net, "local_mean", "global_mean");
    % Connect std branch to embedding layer.
    net = connectLayers(net, "local_std", "global_std");
end

%% Create MAPPO Agents: one actor per UAV sharing the centralized critic.
agentOptions = rlPPOAgentOptions(... 
    'ClipFactor', 0.2, ...
    'EntropyLossWeight', 0.01, ...
    'MiniBatchSize', 256, ...
    'ExperienceHorizon', env.maxSteps, ...
    'DiscountFactor', 0.99);

agents = cell(numUAVs, 1);
for i = 1:numUAVs
    % Create the actor network for UAV i.
    actorNet = createActorNetForAgent(i, numUAVs, stepSize);
    
    % Create the actor representation using a Gaussian policy.
    % Now, the actor outputs a global action vector of size [numUAVs*2,1],
    % with nonzero entries only at the positions corresponding to UAV i.
    actor = rlContinuousGaussianActor(actorNet, globalObsInfo, actorActInfo, ...
        'ObservationInputNames', 'globalState', ...  % Both actor and critic use global observation.
        'ActionMeanOutputNames', 'global_mean', ...
        'ActionStandardDeviationOutputNames', 'global_std');
    actor.UseDevice = "gpu";
    
    % Create the PPO agent for UAV i.
    agents{i} = rlPPOAgent(actor, critic, agentOptions);
end

%% Training Options
trainOpts = rlTrainingOptions(... 
    'MaxEpisodes', 100, ...
    'MaxStepsPerEpisode', env.maxSteps, ...
    'StopTrainingCriteria', 'EpisodeReward', ...
    'StopTrainingValue', Inf, ...
    'Verbose', true, ...
    'Plots', 'training-progress');

%% Train Agents (for simplicity, train each UAV agent sequentially)
for i = 1:numUAVs
    fprintf('Training agent for UAV %d\n', i);
    trainingStats = train(agents{i}, env, trainOpts);
end

%% Save the Trained Agents and Critic
save("MAPPO_UAV_Agents.mat", "agents", "critic");