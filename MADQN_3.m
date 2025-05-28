clc;
clear;
close all;

%% Environment Setup
numUEs = 400;     % Number of UEs
height = 100;
areaSize = 2000;  % Size of area (e.g., 2000x2000)

dataRate = [100*10^6, 15*10^6, 7*10^6];

alpha = 4.88;
beta = 0.43;

[Zcov, Zcap, Rcov, Rcap, ~] = numerology(height,areaSize,areaSize,numUEs, 0.20, dataRate, alpha, beta);

numUAVs = max(Zcov, Zcap); 
Rlimit = min(Rcov, Rcap);

stepSize = 0.3*Rlimit;

% File directory and data loading (update path as needed)
dir = '\\maa1.cc.lut.fi\home\z143799\Documents\Projects\uav-ahp-main\';
UEPositions = load(fullfile(dir, 'Users1.txt'));
demandVector = load(fullfile(dir, 'Demand1.txt'));

% Create the environment
env = MADQN_Environment(numUAVs, numUEs, UEPositions, demandVector, height, Rlimit, areaSize, stepSize);

%% Get Observation and Action Specifications
obsInfo = getObservationInfo(env);  % rlNumericSpec with Dimension [8 1]
actInfo = getActionInfo(env);       % rlFiniteSetSpec for joint actions (1:625)

%% Build Q-Network with Dropout Layers
inputDim = prod(obsInfo.Dimension);
actionInputDim = 1;

% Observation Branch
obsBranch = [
    featureInputLayer(inputDim, "Normalization", "none", "Name", "obs")
    fullyConnectedLayer(128, "Name", "obs_fc1")
    reluLayer("Name", "obs_relu1")
    dropoutLayer(0.2, "Name", "obs_dropout1")
    fullyConnectedLayer(64, "Name", "obs_fc2")
    reluLayer("Name", "obs_relu2")
    dropoutLayer(0.2, "Name", "obs_dropout2")
];

% Action Branch
actBranch = [
    featureInputLayer(actionInputDim, "Normalization", "none", "Name", "act")
    fullyConnectedLayer(16, "Name", "act_fc1")
    reluLayer("Name", "act_relu1")
    dropoutLayer(0.2, "Name", "act_dropout1")
];

% Common Branch
combinedLayer = concatenationLayer(1,2,"Name","concat");
commonLayers = [
    fullyConnectedLayer(64, "Name", "common_fc1")
    reluLayer("Name", "common_relu1")
    dropoutLayer(0.2, "Name", "common_dropout1")
    fullyConnectedLayer(1, "Name", "q_output") % Scalar Q-value output
];

% Assemble Network
lgraph = layerGraph();
lgraph = addLayers(lgraph, obsBranch);
lgraph = addLayers(lgraph, actBranch);
lgraph = addLayers(lgraph, combinedLayer);
lgraph = addLayers(lgraph, commonLayers);

lgraph = connectLayers(lgraph, "obs_dropout2", "concat/in1");
lgraph = connectLayers(lgraph, "act_dropout1", "concat/in2");
lgraph = connectLayers(lgraph, "concat", "common_fc1");

%% Create Q-Value Function Representation
critic = rlQValueFunction(lgraph, obsInfo, actInfo, ...
    "ObservationInputNames", "obs", "ActionInputNames", "act");

%% Move Critic to GPU
critic = moveNetworkToGPU(critic);
agent = rlDQNAgent(critic, rlDQNAgentOptions());

%% Define Training Parameters
maxEpisodes = 1000;
maxSteps = env.maxSteps;
batchSize = 128;
gamma = 0.99;
learningRate = 0.001;

% Start GPU-accelerated training
customTrainMADQN(agent, env, maxEpisodes, maxSteps, batchSize, gamma, learningRate);

%% Save the Trained Agent
save("MADQN_Agent.mat", "agent");

%% Custom GPU Training Function
function customTrainMADQN(agent, env, maxEpisodes, maxSteps, batchSize, gamma, learningRate)
    criticNet = getCritic(agent);
    criticNet = moveNetworkToGPU(criticNet);  % Ensure critic network is moved to GPU
    agent = setCritic(agent, criticNet);

    gradDecay = 0.9;
    squaredGradDecay = 0.99;
    optimizer = adamupdate(learningRate, gradDecay, squaredGradDecay);

    replayBuffer = rlReplayMemory(getObservationInfo(env), getActionInfo(env), 1e6);

    for episode = 1:maxEpisodes
        state = gpuArray(single(reset(env)));
        totalReward = 0;

        for step = 1:maxSteps
            action = getEpsilonGreedyAction(agent, state);

            [nextState, reward, done] = step(env, action);
            nextState = gpuArray(single(nextState));

            experience = {state, action, reward, nextState, done};
            append(replayBuffer, experience);

            if length(replayBuffer) > batchSize
                batch = sample(replayBuffer, batchSize);
                agent = trainStep(agent, batch, gamma, optimizer);
            end

            state = nextState;
            totalReward = totalReward + reward;
            if done, break; end
        end

        fprintf("Episode %d | Total Reward: %.2f\n", episode, totalReward);
    end
end

%% Custom Training Step Function
function agent = trainStep(agent, batch, gamma, optimizer)
    states = gpuArray(batch.Observations);
    actions = gpuArray(batch.Actions);
    rewards = gpuArray(batch.Rewards);
    nextStates = gpuArray(batch.NextObservations);
    dones = gpuArray(batch.IsDone);

    targetQ = rewards + gamma * max(predict(agent.Critic, nextStates), [], 2) .* (1 - dones);
    predictedQ = predict(agent.Critic, states);
    predictedQ = predictedQ(sub2ind(size(predictedQ), (1:length(actions))', actions));

    loss = mean((targetQ - predictedQ).^2);
    gradients = dlgradient(loss, agent.Critic.Learnables);

    [agent.Critic.Learnables, optimizer] = adamupdate(agent.Critic.Learnables, gradients, optimizer);
end

%% Epsilon-Greedy Action Selection
function action = getEpsilonGreedyAction(agent, state)
    persistent epsilon;
    if isempty(epsilon)
        epsilon = 1;  
    end

    if rand < epsilon
        action = randi([1 625]);
    else
        action = getAction(agent, state);
    end

    epsilon = max(0.2, epsilon * 0.999);
end

%% Move Critic to GPU
function critic = moveNetworkToGPU(critic)
    % Move all layers in the critic network to GPU
    layers = critic.LayerGraph.Layers;
    
    % Move the parameters of each layer to the GPU
    for i = 1:numel(layers)
        if isprop(layers{i}, 'Weights')
            layers{i}.Weights = gpuArray(layers{i}.Weights);
        end
        if isprop(layers{i}, 'Bias')
            layers{i}.Bias = gpuArray(layers{i}.Bias);
        end
    end
    
    % Rebuild the layer graph with GPU parameters
    critic.LayerGraph.Layers = layers;
end