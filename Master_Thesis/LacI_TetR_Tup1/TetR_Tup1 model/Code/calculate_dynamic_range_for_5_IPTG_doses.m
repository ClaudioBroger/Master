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
ModelSettings_LacI_TetR_model;

doses = [1 25 50 75 100];

%%Find Baseline expression
dynamic_range_rep1 = zeros(num_rows, width(doses));
for k = 1:num_rows
    for dose = 1:width(doses)
        data_baseline_rep1 = SimFluoValues1_combined_aTc(SimFluoValues1_combined_aTc.draw == k,:);
        total_columns = size(data_baseline_rep1,2);
        data_baseline_rep1(:, (total_columns-2):total_columns) = [];
        
        % Calculate the minimum value
        min_value = min(min(table2array(data_baseline_rep1(:,doses(dose)))));
        baseline_rep1(k,dose) = min_value;
         [row_index, col_index] = find(table2array(data_baseline_rep1(:,doses(dose))) == min_value);
         aTc_baseline_rep1(k,dose) = data.dose(row_index);
         IPTG_baseline_rep1(k,dose) = data_IPTG.dose(doses(dose));


        % Calculate the maximum value
        max_value = max(max(table2array(data_baseline_rep1(:,doses(dose)))));
        max_values_rep1(k,dose) = max_value;
         [row_index, col_index] = find(table2array(data_baseline_rep1(:,doses(dose))) == max_value);
         aTc_max_rep1(k,dose) = data.dose(row_index);
         IPTG_max_rep1(k,dose) = data_IPTG.dose(doses(dose));


        % Calculate the dynamic range
        dynamic_range_rep1(k,dose) = max_value - min_value;

        data_baseline_rep2 = SimFluoValues2_combined_aTc(SimFluoValues2_combined_aTc.draw == k,:);
        total_columns = size(data_baseline_rep2,2);
        data_baseline_rep2(:, (total_columns-2):total_columns) = [];
        
        % Calculate the minimum value
        min_value = min(min(table2array(data_baseline_rep2(:,doses(dose)))));
        baseline_rep2(k,dose) = min_value;
         [row_index, col_index] = find(table2array(data_baseline_rep2(:,doses(dose))) == min_value);
         aTc_baseline_rep2(k,dose) = data.dose(row_index);
         IPTG_baseline_rep2(k,dose) = data_IPTG.dose(doses(dose));


        % Calculate the maximum value
        max_value = max(max(table2array(data_baseline_rep2(:,doses(dose)))));
        max_values_rep2(k,dose) = max_value;
         [row_index, col_index] = find(table2array(data_baseline_rep2(:,doses(dose))) == max_value);
         aTc_max_rep2(k,dose) = data.dose(row_index);
         IPTG_max_rep2(k,dose) = data_IPTG.dose(doses(dose));


        % Calculate the dynamic range
        dynamic_range_rep2(k,dose) = max_value - min_value;

        data_baseline_rep3 = SimFluoValues3_combined_aTc(SimFluoValues3_combined_aTc.draw == k,:);
        total_columns = size(data_baseline_rep3,2);
        data_baseline_rep3(:, (total_columns-2):total_columns) = [];
        
        % Calculate the minimum value
        min_value = min(min(table2array(data_baseline_rep3(:,doses(dose)))));
        baseline_rep3(k,dose) = min_value;
         [row_index, col_index] = find(table2array(data_baseline_rep3(:,doses(dose))) == min_value);
         aTc_baseline_rep3(k,dose) = data.dose(row_index);
         IPTG_baseline_rep3(k,dose) = data_IPTG.dose(doses(dose));


        % Calculate the maximum value
        max_value = max(max(table2array(data_baseline_rep3(:,doses(dose)))));
        max_values_rep3(k,dose) = max_value;
         [row_index, col_index] = find(table2array(data_baseline_rep3(:,doses(dose))) == max_value);
         aTc_max_rep3(k,dose) = data.dose(row_index);
         IPTG_max_rep3(k,dose) = data_IPTG.dose(doses(dose));


        % Calculate the dynamic range
        dynamic_range_rep3(k,dose) = max_value - min_value;
    end
end


for k = 1:num_rows
    figure(k);
    marker_size = 40; % Adjust the marker size as needed
    h = scatter3(IPTG_baseline_rep1(k,:), aTc_baseline_rep1(k,:), baseline_rep1(k,:), marker_size, baseline_rep1(k,:), 'filled');
    xlabel('IPTG Doses');
    ylabel('aTc Doses');
    zlabel('Baseline Rep1');
    title('Baseline Rep 1 of IPTG and aTc Doses');
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Baseline Rep 1'); % Label the colorbar
    grid on;
end

figure;
plot(IPTG_baseline_rep1, dynamic_range_rep1)

for k = 1:num_rows
    figure(k);
    marker_size = 40; % Adjust the marker size as needed
    h = scatter3(IPTG_baseline_rep2(k,:), aTc_baseline_rep2(k,:), baseline_rep2(k,:), marker_size, baseline_rep2(k,:), 'filled');
    xlabel('IPTG Doses');
    ylabel('aTc Doses');
    zlabel('Baseline Rep 2');
    title('Baseline Rep 2 of IPTG and aTc Doses');
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Baseline Rep 2'); % Label the colorbar
    grid on;
end

for k = 1:num_rows
    figure(k);
    marker_size = 40; % Adjust the marker size as needed
    h = scatter3(IPTG_baseline_rep3(k,:), aTc_baseline_rep3(k,:), baseline_rep3(k,:), marker_size, baseline_rep2(k,:), 'filled');
    xlabel('IPTG Doses');
    ylabel('aTc Doses');
    zlabel('Dynamic Range');
    title('Dynamic Range of IPTG and aTc Doses');
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Dynamic Range'); % Label the colorbar
    grid on;
end

for k = 1:num_rows
    figure(k);
    marker_size = 40; % Adjust the marker size as needed
    h = scatter3(IPTG_max_rep1(k,:), aTc_max_rep1(k,:), max_values_rep1(k,:), marker_size, max_values_rep1(k,:), 'filled');
    xlabel('IPTG Doses');
    ylabel('aTc Doses');
    zlabel('Dynamic Range');
    title('Dynamic Range of IPTG and aTc Doses');
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Dynamic Range'); % Label the colorbar
    grid on;
end

for k = 1:num_rows
    figure(k);
    marker_size = 40; % Adjust the marker size as needed
    h = scatter3(IPTG_max_rep2(k,:), aTc_max_rep2(k,:), max_values_rep2(k,:), marker_size, max_values_rep2(k,:), 'filled');
    xlabel('IPTG Doses');
    ylabel('aTc Doses');
    zlabel('Dynamic Range');
    title('Dynamic Range of IPTG and aTc Doses');
    colormap(jet); % Or any other colormap of your choice
    c = colorbar; % Display the colorbar to show the mapping of dynamic range to colors
    ylabel(c, 'Dynamic Range'); % Label the colorbar
    grid on;
end


