load('Time_dynamics_fluo_over_IPTG_rep1.mat')
load('Time_dynamics_fluo_over_IPTG_rep2.mat')
load('Time_dynamics_fluo_over_IPTG_rep3.mat')
load('Time_IPTG_rep1.mat')
load('Time_IPTG_rep2.mat')
load('Time_IPTG_rep3.mat')

sbioreset;
simple_model_LacI_TetR;
sbiosaveproject('simple_model');
model = sbioloadproject('simple_model');

sbioaccelerate(model.mw_sbmod1)

%Modelsettings
ModelSettings_LacI_TetR;

%Time dynamics for growth rate over IPTG dose Rep 1
C = linspecer(height(data.dose));
for k = 1:width(SimFluoValues_time1_combined)
    figure(k);
     for a = 1:height(data.dose)


        hold on;

        rounded_dose = round(data.dose(a) * 100) / 100;

        plot(time_IPTG1_combined{a,k},SimFluoValues_time1_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 50)
        ylabel('mean Fluorescence','FontSize', 50)
        title('Rep 1 mean Fluorescence IPTG time dynamics', 'FontSize',60)

    
    end
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size         

    colors = C; 
    colormap(colors);

% Create colorbar
    cb = colorbar;

% Set colorbar ticks and labels
    cb.Ticks = [0 0.5 1]; % positions of the midpoints of the color segments
    cb.TickLabels = {data.dose(1), data.dose(500), data.dose(1000)};
    cb.FontSize = 25;
    ylabel(cb, 'IPTG dose','FontSize', 30)


end
close all

%Time dynamics for growth rate over IPTG dose Rep 2
C = linspecer(height(data.dose));
for k = 1:width(SimFluoValues_time2_combined)
    figure(k);
     for a = 1:height(data.dose)


        hold on;

        rounded_dose = round(data.dose(a) * 100) / 100;

        plot(time_IPTG2_combined{a,k},SimFluoValues_time2_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 50)
        ylabel('mean Fluorescence','FontSize', 50)
        title('Rep 2 mean Fluorescence IPTG time dynamics', 'FontSize',60)

    
    end
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size         

    colors = C; 
    colormap(colors);

% Create colorbar
    cb = colorbar;

% Set colorbar ticks and labels
    cb.Ticks = [0 0.5 1]; % positions of the midpoints of the color segments
    cb.TickLabels = {data.dose(1), data.dose(500), data.dose(1000)};
    cb.FontSize = 25;
    ylabel(cb, 'IPTG dose','FontSize', 30)



end
close all

%Time dynamics for growth rate over IPTG dose Rep 3
C = linspecer(height(data.dose));
for k = 1:width(SimFluoValues_time3_combined)
    figure(k);
     for a = 1:height(data.dose)


        hold on;

        rounded_dose = round(data.dose(a) * 100) / 100;

        plot(time_IPTG3_combined{a,k},SimFluoValues_time3_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 50)
        ylabel('mean Fluorescence','FontSize', 50)
        title('Rep 3 mean Fluorescence IPTG time dynamics', 'FontSize',60)

    
    end
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size         

    colors = C; 
    colormap(colors);

% Create colorbar
    cb = colorbar;

% Set colorbar ticks and labels
    cb.Ticks = [0 0.5 1]; % positions of the midpoints of the color segments
    cb.TickLabels = {data.dose(1), data.dose(500), data.dose(1000)};
    cb.FontSize = 25;
    ylabel(cb, 'IPTG dose','FontSize', 30)


end
close all





