classdef LocalStateExtractionLayer < nnet.layer.Layer
    properties
        AgentIndex
    end
    methods
        function layer = LocalStateExtractionLayer(agentIndex, name)
            % Set layer properties.
            layer.AgentIndex = agentIndex;
            layer.Name = name;
            layer.Description = "Extracts UAV " + num2str(agentIndex) + " local state from global state";
        end
        
        function Z = predict(layer, X)
            % X is expected to be of size [numUAVs*2, N] where N is the mini-batch size.
            startIdx = 2 * layer.AgentIndex - 1;
            endIdx = 2 * layer.AgentIndex;
            Z = X(startIdx:endIdx, :);
        end
    end
end
