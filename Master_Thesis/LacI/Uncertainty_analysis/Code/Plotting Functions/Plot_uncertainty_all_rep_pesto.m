dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
dataPos = strcat(dataPath, "data.mat");
data1 = load(dataPos);
data1.dose = [linspace(min(data1.dose), 10^3, 250) linspace(1+10^3, 10^7, 10000) linspace(1+10^7, max(data1.dose), 100)];
data1.dose = data1.dose';

dataPos = strcat(dataPath, "data.mat");
data2 = load(dataPos);
data2.dose = [linspace(min(data1.dose), 10^5, 100) linspace(1+10^5, 10^7, 1500) linspace(1+10^7, max(data1.dose), 100)];
data2.dose = data2.dose';

dataPos = strcat(dataPath, "data.mat");
data3 = load(dataPos);
data3.dose = [linspace(min(data1.dose), 10^4, 100) linspace(1+10^4, 10^9, 10000) linspace(1+10^9, 10^10, 100)];
data3.dose = data3.dose';








figure(1)
hold on

% Let's create empty arrays to store the plot handles
p1 = [];
p2 = [];
p3 = [];

for draw = 1:30
    % For the first iteration, save the plot handles
    if draw == 1
        p1 = plot(log10(data1.dose),SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'red');
        p2 = plot(log10(data2.dose),SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'blue');
        p3 = plot(log10(data3.dose),SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'green');
    else % For other iterations, use 'handlevisibility' to 'off'
        plot(log10(data1.dose),SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'red', 'HandleVisibility','off');
        plot(log10(data2.dose),SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'blue', 'HandleVisibility','off');
        plot(log10(data3.dose),SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'green', 'HandleVisibility','off');
    end
end

% Create the legend after the loop
legend([p1 p2 p3],{'Repression coefficient 1', 'Repression coefficient 2', 'Repression coefficient 3'},'Location', 'northeastoutside','FontSize',15);
xlabel('log IPTG (nM)', 'FontSize', 30)
ylabel('mean Fluorescence','FontSize', 30)
title('Uncertainty all Repression coefficients - PESTO', 'FontSize',40)
      