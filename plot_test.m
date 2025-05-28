figure
hold on

% === Settings ===
numScenarios = 3;
numPlans = 3;
numMethods = 4;

scenarioLabels = {'Scenario 1', 'Scenario 2', 'Scenario 3'};
planLabels = {'SINR', 'TP', 'Coverage'};
methodLabels = {'UAV-AHP', 'PSO', 'CS', 'NSGA-II'};
methodColors = [0, 114, 178; 240, 228, 66; 213, 94, 0; 86, 180, 233]/255;

% === Simulated Data (replace with your real values) ===
Data = rand(numScenarios, numPlans, numMethods)*0.1 + 0.85;        % Mean values
Errors = rand(numScenarios, numPlans, numMethods)*0.01 + 0.005;    % Error bars

% === Layout Settings ===
methodSpacing = 1;
planSpacing = 5;        % Space between plans
scenarioSpacing = 2;    % Space between scenarios

% === Generate y positions ===
y = [];
yLabels = {};
barIdx = 0;
yGroupCenters = zeros(numScenarios, numPlans);  % To store center positions of each plan group

for s = 1:numScenarios
    for p = 1:numPlans
        baseY = (s-1)*(numPlans*planSpacing + scenarioSpacing) + (p-1)*planSpacing;
        for m = 1:numMethods
            barIdx = barIdx + 1;
            y(barIdx) = baseY + m*methodSpacing;
            yLabels{barIdx} = sprintf('%s - %s - %s', scenarioLabels{s}, planLabels{p}, methodLabels{m});
        end
        % Save center of this plan's group (for labeling later)
        yGroupCenters(s, p) = baseY + (numMethods + 1)/2 * methodSpacing;
    end
end

% === Flatten data for plotting ===
flatData = reshape(permute(Data, [3 2 1]), [], 1);
flatErrors = reshape(permute(Errors, [3 2 1]), [], 1);

% === Plot bars and error bars ===
for i = 1:length(flatData)
    methodIdx = mod(i-1, numMethods) + 1;
    barh(y(i), flatData(i), 0.8, 'FaceColor', methodColors(methodIdx,:), 'EdgeColor', 'none');
    errorbar(flatData(i), y(i), flatErrors(i), 'horizontal', '.', 'Color', 'k', 'LineWidth', 1);
end

% === Add dashed lines between scenarios ===
maxY = max(y) + 3;
for s = 1:numScenarios-1
    splitPos = s*(numPlans*planSpacing + scenarioSpacing) - scenarioSpacing/2;
    xLims = xlim;
    plot(xLims, [splitPos splitPos], 'k--', 'LineWidth', 1);
end

% === Label plans on Y axis ===
% Replace YTick labels with plan names at each plan center
yticks(reshape(yGroupCenters', [], 1));
yticklabels(repmat(planLabels, 1, numScenarios));

% === X-axis settings ===
xlabel('Mean Coverage of Users')
xlim([0.6 1.02])
ylim([0 maxY])
set(gca, 'YDir','reverse') % So Scenario 1 is on top (optional)

% === Add Scenario Labels ===
for s = 1:numScenarios
    % Compute vertical center of this scenario's full group
    groupTop = (s-1)*(numPlans*planSpacing + scenarioSpacing);
    groupBottom = groupTop + numPlans*planSpacing;
    groupCenterY = (groupTop + groupBottom) / 2;

    text(0.585, groupCenterY, scenarioLabels{s}, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 10)
end

% % === Add Plan Labels (SINR, TP, Coverage) to the right of the ticks ===
% for s = 1:numScenarios
%     for p = 1:numPlans
%         planY = yGroupCenters(s, p);
%         text(1.025, planY, planLabels{p}, ...
%             'HorizontalAlignment', 'left', ...
%             'VerticalAlignment', 'middle', ...
%             'FontAngle', 'italic', 'FontSize', 9)
%     end
% end

% === Legend and Title ===
legend(methodLabels, 'Location', 'southoutside', 'Orientation', 'horizontal')
title('Coverage by Scenario, Plan, and Method')
box on
hold off
