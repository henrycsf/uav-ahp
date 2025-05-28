classdef MAPPO_Environment < rl.env.MATLABEnvironment
    properties
        numUAVs
        numUEs
        height
        areaSizeX      % Maximum x-value of the area
        areaSizeY      % Maximum y-value of the area
        stepSize       % 2-element vector: [stepSizeX, stepSizeY]
        UAVPositions
        UEPositions
        initialPositions
        initialCoverage
        initialSINR
        initialTP
        initialDemand
        cellRadius
        demandVector
        maxSteps = 10;
        SINRThreshold = -10 % SINR minimum threshold in dB
    end
    
    properties (Access = private)
        currentStep = 0;
    end

    methods
        function obj = MAPPO_Environment(numUAVs, numUEs, UEPositions, demandVector, height, areaSizeX, areaSizeY, stepSize)
            % Define observation dimension: each UAV has (x,y) so total dimension = numUAVs*2.
            obsDim = numUAVs * 2;
            % Create upper limit vector by repeating [areaSizeX; areaSizeY] for each UAV.
            upperLimit = repmat([areaSizeX(1); areaSizeY(1)], numUAVs, 1);
            observationInfo = rlNumericSpec([obsDim 1], 'LowerLimit', zeros(obsDim,1), 'UpperLimit', upperLimit);
            observationInfo.Name = 'state';
            
            % For continuous actions, each UAV chooses a 2-D movement vector.
            actionDim = numUAVs * 2;
            % Allowable action: each UAV's dx is in [-stepSizeX, stepSizeX] and dy in [-stepSizeY, stepSizeY]
            lowerLimit = -repmat(stepSize(:), numUAVs, 1);
            upperLimit = repmat(stepSize(:), numUAVs, 1);
            actionInfo = rlNumericSpec([actionDim 1], 'LowerLimit', lowerLimit, 'UpperLimit', upperLimit);
            actionInfo.Name = 'action';
            
            % Call superclass constructor.
            obj@rl.env.MATLABEnvironment(observationInfo, actionInfo);

            load("numerology_festival.mat", "CaseSelect", "Pt", "DataRate", "Rcap", "Rcov", "Zcap", "Zcov", "nmi","alpha","beta");

            plots = false;

            [UAV_info, Best, ~, ~] = drone_positioning_AHP(CaseSelect, Pt, numUEs, areaSizeX, ...
                areaSizeY, height, DataRate, UEPositions, demandVector, Zcov, Zcap, Rcov, Rcap, nmi, alpha, beta, plots);

            for i = 1:numUAVs
                UAVPositions(i,:) = [UAV_info(Best(:,2,i),1,i), UAV_info(Best(:,2,i),2,i)];
            end

            cellRadius = min(Rcov, Rcap)*10^3;
            
            % Assign properties.
            obj.numUAVs = numUAVs;
            obj.numUEs = numUEs;
            obj.height = height;
            obj.areaSizeX = areaSizeX;
            obj.areaSizeY = areaSizeY;
            obj.stepSize = stepSize;  % e.g., [stepSizeX, stepSizeY]
            obj.demandVector = demandVector;
            obj.cellRadius = cellRadius;
            % Initialize UAV positions randomly within the rectangular area.
            obj.UAVPositions = UAVPositions;
            obj.initialPositions = UAVPositions;
            % Assume UEPositions and demandVector are provided externally.
            obj.UEPositions = UEPositions;
        end
        
        function state = reset(obj)
            % Optionally, you could use clustering to set initial UAV positions.
            % Here we reinitialize UAV positions to the initial positions.

            load("numerology_festival.mat", "CaseSelect", "Pt", "DataRate", "Rcap", "Rcov", "Zcap", "Zcov", "nmi","alpha","beta");

            plots = false;

            [UAV_info, Best, ~, ~] = drone_positioning_AHP(CaseSelect, Pt, obj.numUEs, obj.areaSizeX, ...
                obj.areaSizeY, obj.height, DataRate, obj.UEPositions, obj.demandVector, Zcov, Zcap, Rcov, Rcap, nmi, alpha, beta, plots);
            
            for i = 1:obj.numUAVs
                obj.UAVPositions(i,:) = [UAV_info(Best(:,2,i),1,i), UAV_info(Best(:,2,i),2,i)];
            end

            % Optionally, add noise to UE positions to simulate dynamics.
            obj.UEPositions = obj.UEPositions + (5*randn(obj.numUEs,2));
            obj.currentStep = 0;
            state = obj.UAVPositions(:); % Return as a column vector.
        end
        
        function [nextState, reward, done] = step(obj, actions)
            % Actions is a vector of size [numUAVs*2, 1] representing [dx;dy] for each UAV.
            delta = reshape(actions, obj.numUAVs, 2);
            
            % Update UAV positions:
            newPositions = obj.UAVPositions + delta;
            
            % Reflective clamping: If a UAV overshoots, reflect it back inside.
            for i = 1:obj.numUAVs
                % For x-coordinate
                if newPositions(i,1) < 0
                    newPositions(i,1) = -newPositions(i,1);
                elseif newPositions(i,1) > obj.areaSizeX
                    newPositions(i,1) = obj.areaSizeX - (newPositions(i,1)-obj.areaSizeX);
                end
                % For y-coordinate
                if newPositions(i,2) < 0
                    newPositions(i,2) = -newPositions(i,2);
                elseif newPositions(i,2) > obj.areaSizeY
                    newPositions(i,2) = obj.areaSizeY - (newPositions(i,2)-obj.areaSizeY);
                end
            end
            obj.UAVPositions = newPositions;
            
            % Compute performance metrics:
            [coverage, SINR, throughput, demand] = obj.evaluateCoverage();
            
            % Normalize performance metrics (scaling factors may need tuning)
            normalCoverage = coverage / obj.initialCoverage;
            normalSINR = abs((SINR) / (obj.initialSINR));   % Here, SINR is normalized from -10 to 0 dB.
            normalTP = throughput / obj.initialTP;
            normalDemand = demand / obj.initialDemand;
            
            % Increase punishment for UAVs staying at or near boundaries.
            punishment = 0;
            edgeThresholdX = 0.01 * obj.areaSizeX;  % 1% of area width.
            edgeThresholdY = 0.01 * obj.areaSizeY;  % 1% of area height.
            for i = 1:obj.numUAVs
                if obj.UAVPositions(i,1) <= edgeThresholdX || obj.UAVPositions(i,1) >= (obj.areaSizeX - edgeThresholdX)
                    punishment = punishment + 0.5;
                end
                if obj.UAVPositions(i,2) <= edgeThresholdY || obj.UAVPositions(i,2) >= (obj.areaSizeY - edgeThresholdY)
                    punishment = punishment + 0.5;
                end
            end
            
            % Adjust reward formulation as needed. Here we subtract the increased punishment.
            reward = (10*normalCoverage + 10*normalSINR + 5*normalTP + 2*normalDemand - 10*punishment);
            
            nextState = obj.UAVPositions(:);
            obj.currentStep = obj.currentStep + 1;
            done = (obj.currentStep >= obj.maxSteps);
        end

        
        function [coverage, SINR, throughput, demand] = evaluateCoverage(obj)
            % This function computes network performance metrics based on current UAV positions.
            f = 3.5 * 10^9;
            velc = 299792458;
            ZetaLOS = 1;
            ZetaNLOS = 20;
            BW = [50*10^6, 20*10^6, 10*10^6];
            q = -174 + 10*log10(BW);
            Gt = 3;
            Gr = 0;
            h = obj.height;
            Pt = 35;
            alpha = 4.88;
            beta = 0.43;
            BW = 6*10^6;
            
            SINR = 0;
            coverage = 0;
            throughput = 0;
            demand = 0;
            
            for ue = 1:obj.numUEs
                R = vecnorm(obj.UAVPositions - obj.UEPositions(ue,:), 2, 2);
                demandDensity = obj.demandVector(ue) ./ (R.^2);
                theta = atan(h./R);
                Z = (alpha*exp(-beta*((180/pi).*theta - alpha)));
                PL = 20*log10((4*pi*f*R./velc)) + ((ZetaLOS+Z.*ZetaNLOS)./(1+Z));
                signalStrength = Pt - PL + Gt + Gr;
                linearSSRI = (10.^((signalStrength-30)./10));
                interference = sum(linearSSRI) - max(linearSSRI);
                maxLinearSINR = max(linearSSRI) / ((10^((q(1)-30)/10)) + interference);
                maxLinearSNR = max(linearSSRI) / (10^((q(1)-30)/10));
                maxSINR = 10 * log10(maxLinearSINR);
                maxTP = (10^-6)*BW*log2(1+maxLinearSNR);
                if maxSINR > obj.SINRThreshold && min(vecnorm(obj.UAVPositions - obj.UEPositions(ue,:), 2, 2)) <= obj.cellRadius
                    coverage = coverage + 1;
                    SINR = SINR + maxSINR;
                    throughput = throughput + maxTP;
                    demand = demand + max(demandDensity);

                else
                    SINR = + obj.SINRThreshold;
                    demand = demand + max(demandDensity);
                end
            end
            coverage = coverage / obj.numUEs;
            SINR = SINR / obj.numUEs;
            throughput = throughput / obj.numUEs;
            demand = demand / obj.numUEs;

            if obj.currentStep == 0
                obj.initialCoverage = coverage;
                obj.initialSINR = SINR;
                obj.initialTP = throughput;
                obj.initialDemand = demand;
            end
        end
    end
end
