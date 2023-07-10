figure(1)
hold on

% Let's create empty arrays to store the plot handles
p1 = [];
p2 = [];
p3 = [];

for draw = 1:30
    % For the first iteration, save the plot handles
    if draw == 1
        p1 = plot(log10(data.dose),SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'red');
        p2 = plot(log10(data.dose),SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'blue');
        p3 = plot(log10(data.dose),SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'green');
    else % For other iterations, use 'handlevisibility' to 'off'
        plot(log10(data.dose),SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'red', 'HandleVisibility','off');
        plot(log10(data.dose),SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'blue', 'HandleVisibility','off');
        plot(log10(data.dose),SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == draw),'-', 'LineWidth', 2, 'Color', 'green', 'HandleVisibility','off');
    end
end

% Create the legend after the loop
legend([p1 p2 p3],{'Repression coefficient 1', 'Repression coefficient 2', 'Repression coefficient 3'},'Location', 'northeastoutside','FontSize',30);
xlabel('log IPTG (nM)', 'FontSize', 60)
ylabel('mean Fluorescence','FontSize', 60)
title('LacI TetRTup1 Uncertainty all Repression coefficients - no aTc', 'FontSize',70)
ax = gca; % Get the current axes
ax.XAxis.FontSize = 40; % X axis font size
ax.YAxis.FontSize = 40; % Y axis font size   

set(figure(1), 'Position', get(0, 'Screensize'));
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'LacI_TetRTup1_uncertainty_all_rep_no_aTc', '.jpg']);
saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'LacI_TetRTup1_uncertainty_all_rep_no_aTc', '.fig']);

close all