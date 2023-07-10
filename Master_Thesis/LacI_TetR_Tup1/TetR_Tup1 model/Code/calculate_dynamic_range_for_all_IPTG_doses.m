%%Calculate dynamic range of selected simulations

load('05-May-2023SimFluoValues1_combined_LacI_TetRTup1.mat');
load('05-May-2023SimFluoValues2_combined_LacI_TetRTup1.mat');
load('05-May-2023SimFluoValues3_combined_LacI_TetRTup1.mat');

load("05-May-2023_rand_parameter_selection2_LacI_TetRTup1.mat")
rand_parameter1 = load('05-May-2023_rand_parameter_selection_LacI_TetRTup1.mat');
rand_parameter1 = rand_parameter1.rand_parameter;
rand_parameter1(1,:) = [];
rand_parameter = [rand_parameter; rand_parameter1];

num_rows = height(rand_parameter);

%Reset, reload and compile model
sbioreset;
LacI_TetRTup1_model_with_aTc;
sbiosaveproject('LacI_TetRTup1_model');
model = sbioloadproject('LacI_TetRTup1_model');


sbioaccelerate(model.mw_sbmod1)
%Modelsettings
ModelSettings_LacI_TetR_model_aTc;

%%Find Baseline expression
dynamic_range_rep1 = zeros(num_rows, height(data_IPTG.dose));
dynamic_range_rep2 = zeros(num_rows, height(data_IPTG.dose));
dynamic_range_rep3 = zeros(num_rows, height(data_IPTG.dose));
for k = 1:num_rows
    for dose = 1:height(data_IPTG.dose)
        data_baseline_rep1 = SimFluoValues1_combined_aTc(SimFluoValues1_combined_aTc.draw == k,:);
        total_columns = size(data_baseline_rep1,2);
        data_baseline_rep1(:, (total_columns-2):total_columns) = [];
        
        % Calculate the minimum value
        min_value = min(min(table2array(data_baseline_rep1(:,dose))));
        baseline_rep1(k,dose) = min_value;
        [row_index, col_index] = find(table2array(data_baseline_rep1(:,dose)) == min_value);
        %aTc_baseline_rep1(k,dose) = data.dose(row_index(1));
        IPTG_baseline_rep1(k,dose) = data_IPTG.dose(col_index(1));


        % Calculate the maximum value
        max_value = max(max(table2array(data_baseline_rep1(:,dose))));
        max_values_rep1(k,dose) = max_value;
        [row_index, col_index] = find(table2array(data_baseline_rep1(:,dose)) == max_value);
        %aTc_max_rep1(k,dose) = data.dose(row_index(1));
        IPTG_max_rep1(k,dose) = data_IPTG.dose(col_index(1));


        % Calculate the dynamic range
         dynamic_range_rep1(k,dose) = max_value / min_value;


        data_baseline_rep2 = SimFluoValues2_combined_aTc(SimFluoValues2_combined_aTc.draw == k,:);
        total_columns = size(data_baseline_rep2,2);
        data_baseline_rep2(:, (total_columns-2):total_columns) = [];
        
        % Calculate the minimum value
        min_value = min(min(table2array(data_baseline_rep2(:,dose))));
        baseline_rep2(k,dose) = min_value;
        [row_index, col_index] = find(table2array(data_baseline_rep2(:,dose)) == min_value);
        %aTc_baseline_rep2(k,dose) = data.dose(row_index(1));
        IPTG_baseline_rep2(k,dose) = data_IPTG.dose(col_index(1));


        % Calculate the maximum value
        max_value = max(max(table2array(data_baseline_rep2(:,dose))));
        max_values_rep2(k,dose) = max_value;
        [row_index, col_index] = find(table2array(data_baseline_rep2(:,dose)) == max_value);
        %aTc_baseline_rep2(k,dose) = data.dose(row_index(1));
        IPTG_baseline_rep2(k,dose) = data_IPTG.dose(col_index(1));


        % Calculate the dynamic range
         dynamic_range_rep2(k,dose) = max_value / min_value;


        data_baseline_rep3 = SimFluoValues3_combined_aTc(SimFluoValues3_combined_aTc.draw == k,:);
        total_columns = size(data_baseline_rep3,2);
        data_baseline_rep3(:, (total_columns-2):total_columns) = [];
        
        % Calculate the minimum value
        min_value = min(min(table2array(data_baseline_rep3(:,dose))));
        baseline_rep3(k,dose) = min_value;
        [row_index, col_index] = find(table2array(data_baseline_rep3(:,dose)) == min_value);
        %aTc_baseline_rep3(k,dose) = data.dose(row_index(1));
        IPTG_baseline_rep3(k,dose) = data_IPTG.dose(col_index(1));


        % Calculate the maximum value
        max_value = max(max(table2array(data_baseline_rep3(:,dose))));
        max_values_rep3(k,dose) = max_value;
        [row_index, col_index] = find(table2array(data_baseline_rep3(:,dose)) == max_value);
        %aTc_baseline_rep3(k,dose) = data.dose(row_index(1));
        IPTG_baseline_rep3(k,dose) = data_IPTG.dose(col_index(1));


        % Calculate the dynamic range
         dynamic_range_rep3(k,dose) = max_value / min_value;

    end
end




for k = 1:num_rows
    figure(k);
    h = scatter3(log10(data_IPTG.dose), data.dose, dynamic_range_rep1(k,:), 40, dynamic_range_rep1(k,:), 'filled');
    xlabel('IPTG Doses', 'FontSize', 18);
    ylabel('aTc Doses','FontSize',18);
    zlabel('Dynamic Range Rep 1', 'FontSize',18);
    title(strcat('Dynamic Range Rep 1 - parameter set: ', num2str(k)), "FontSize",20);
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Dynamic Range Rep 1', 'FontSize', 14); % Label the colorbar
    view(-10, 17);
end

for k = 1:num_rows
    figure(k);
    h = scatter3(log10(data_IPTG.dose), data.dose, dynamic_range_rep2(k,:), 40, dynamic_range_rep2(k,:), 'filled');
    xlabel('IPTG Doses', 'FontSize', 18);
    ylabel('aTc Doses','FontSize',18);
    zlabel('Dynamic Range Rep 2', 'FontSize',18);
    title(strcat('Dynamic Range Rep2 - parameter set: ', num2str(k)), "FontSize",20);
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Dynamic Range Rep 2', 'FontSize', 14); % Label the colorbar
    view(-10, 17);
end

for k = 1:num_rows
    figure(k);
    h = scatter3(log10(data_IPTG.dose), data.dose, dynamic_range_rep2(k,:), 40, dynamic_range_rep2(k,:), 'filled');
    xlabel('IPTG Doses', 'FontSize', 18);
    ylabel('aTc Doses','FontSize',18);
    zlabel('Dynamic Range Rep 3', 'FontSize',18);
    title(strcat('Dynamic Range Rep 3 - parameter set: ', num2str(k)), "FontSize",20);
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Dynamic Range Rep 3', 'FontSize', 14); % Label the colorbar
    view(-10, 17);
end

for k = 1:num_rows
    figure(1)
    hold on;
    plot(log10(data_IPTG.dose),dynamic_range_rep1(k,:),'-', 'LineWidth', 2);
    %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
    xlabel('log IPTG (nM)', 'FontSize', 30)
    ylabel('Dynamic range','FontSize', 30)
    title('Dynamic range LacI TetR-Tup1 repression coefficient 1', 'FontSize',40)
end



for k = 1:num_rows
    figure(2)
    hold on;
    plot(log10(data_IPTG.dose),dynamic_range_rep2(k,:),'-', 'LineWidth',2);
    %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
    xlabel('log IPTG (nM)', 'FontSize', 30)
    ylabel('Dynamic range','FontSize', 30)
    title('Dynamic range LacI TetR-Tup1 repression coefficient 2', 'FontSize',40)
end



for k = 1:num_rows
    figure(3)
    hold on;
    plot(log10(data_IPTG.dose),dynamic_range_rep3(k,:),'-', 'LineWidth', 2);
    %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
    xlabel('log IPTG (nM)', 'FontSize', 30)
    ylabel('Dynamic range','FontSize', 30)
    title('Dynamic range LacI TetR-Tup1 repression coefficient 3', 'FontSize',40)
end

figure(4)
hold on;
plot(log10(data_IPTG.dose), dynamic_range_rep1(28,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 1', 'Color', 'red');
plot(log10(data_IPTG.dose), dynamic_range_rep2(28,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 2', 'Color', 'blue');
plot(log10(data_IPTG.dose), dynamic_range_rep3(28,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 3', 'Color', 'green');
xlabel('log IPTG (nM)', 'FontSize', 50)
ylabel('Dynamic range','FontSize', 50)
title('Dynamic range LacI TetR-Tup1 all Rep', 'FontSize',60)
legend("show", 'Location', 'northeastoutside', 'FontSize', 25);
ax = gca; % Get the current axes
ax.XAxis.FontSize = 25; % X axis font size
ax.YAxis.FontSize = 25; % Y axis font size

figure(5)
hold on;
plot(log10(data_IPTG.dose), dynamic_range_rep1(29,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 1', 'Color', 'red');
plot(log10(data_IPTG.dose), dynamic_range_rep2(29,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 2', 'Color', 'blue');
plot(log10(data_IPTG.dose), dynamic_range_rep3(29,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 3', 'Color', 'green');
xlabel('log IPTG (nM)', 'FontSize', 50)
ylabel('Dynamic range','FontSize', 50)
title('Dynamic range LacI TetR-Tup1 all Rep', 'FontSize',60)
legend("show", 'Location', 'northeastoutside', 'FontSize', 25);
ax = gca; % Get the current axes
ax.XAxis.FontSize = 25; % X axis font size
ax.YAxis.FontSize = 25; % Y axis font size

figure(6)
hold on;
plot(log10(data_IPTG.dose), dynamic_range_rep1(30,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 1', 'Color', 'red');
plot(log10(data_IPTG.dose), dynamic_range_rep2(30,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 2', 'Color', 'blue');
plot(log10(data_IPTG.dose), dynamic_range_rep3(30,:),'-', 'LineWidth', 2, 'DisplayName', 'Repression coefficient 3', 'Color', 'green');
xlabel('log IPTG (nM)', 'FontSize', 50)
ylabel('Dynamic range','FontSize', 50)
title('Dynamic range LacI TetR-Tup1 all Rep', 'FontSize',60)
legend("show", 'Location', 'northeastoutside', 'FontSize', 25);
ax = gca; % Get the current axes
ax.XAxis.FontSize = 25; % X axis font size
ax.YAxis.FontSize = 25; % Y axis font size