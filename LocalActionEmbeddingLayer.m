classdef LocalActionEmbeddingLayer < nnet.layer.Layer
    properties
        AgentIndex   % The index of the UAV (agent) for which to embed the action
        NumAgents    % Total number of agents (UAVs)
    end
    
    methods
        function layer = LocalActionEmbeddingLayer(agentIndex, numAgents, name)
            % layer = LocalActionEmbeddingLayer(agentIndex, numAgents, name)
            % agentIndex: index of the UAV (1, 2, ..., NumAgents)
            % numAgents: total number of agents
            % name: layer name
            layer.AgentIndex = agentIndex;
            layer.NumAgents = numAgents;
            layer.Name = name;
            layer.Description = "Embeds local action for agent " + num2str(agentIndex) + " into global action vector";
        end
        
        function Z = predict(layer, X)
            % X is a dlarray of size [2, N] (local action for one agent)
            % We output a dlarray Z of size [numAgents*2, N] with zeros everywhere
            % except rows corresponding to the agent's local action.
            sz = size(X);
            N = sz(2); % mini-batch size
            globalDim = layer.NumAgents * 2;
            % Initialize output with zeros of same type as X.
            Z = zeros(globalDim, N, 'like', X);
            % Determine the rows corresponding to the given agent.
            startIdx = 2 * layer.AgentIndex - 1;
            endIdx = startIdx + 1;
            % Insert the local action into the correct positions.
            Z(startIdx:endIdx, :) = X;
        end
    end
end