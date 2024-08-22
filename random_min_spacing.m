function [Points] = random_min_spacing(maxPoints, dist, vis)

% clc;    % Clear the command window.
% close all;  % Close all figures (except those of imtool.)
% clear;  % Erase all existing variables. Or clearvars if you want.
% workspace;  % Make sure the workspace panel is showing.
% format short g;
% format compact;
% vis = true;

fontSize = 20;

% Define how many points you want to place.
% maxPoints = 300;
% Define counter for how many actually get placed.
numPoints = 0;
% Define fail safe, how many iterations you want to keep trying before quitting.
maxIterations = maxPoints * 10;
loopCounter = 1;
% Define how close points can come to other points before being rejected.
% minClosestDistance = 0.02;
% Declare arrays to hold the x and y coordinate values.
x = nan(1, numPoints);
y = nan(1, numPoints);
while numPoints < maxPoints && loopCounter < maxIterations
	% Get a random point.
	xPossible = rand();
	yPossible = rand();
	if numPoints == 0
		% First point automatically is valid.
		numPoints = numPoints + 1;
		x(numPoints) = xPossible;
		y(numPoints) = yPossible;
		continue;
	end
	% Find distances between this point and all others.
	distances = sqrt((x-xPossible) .^ 2 + (y - yPossible) .^ 2);
	if min(distances) >= dist
		% It's far enough away from all the other points, or it's the first point.
		% Include it in our list of acceptable points.
		numPoints = numPoints + 1;
		x(numPoints) = xPossible;
		y(numPoints) = yPossible;
	end
	% Increment the loop counter.
	loopCounter = loopCounter + 1;
end
% Crop away points we were not able to get.
x = x(1:numPoints); % Crop to valid points.
y = y(1:numPoints); % Crop to valid points.

Points = transpose([x; y]);

if vis==true
    
% Plot the points we were able to get.
plot(x, y, 'b.', 'MarkerSize', 15);
grid on;
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
title('Demo of Random Points with Min Distance', 'FontSize', fontSize);

% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 

% Tell user what we were a ble to achieve.
message = sprintf('Found %d points in %d tries', numPoints, loopCounter);
helpdlg(message);

end
