load('27-May-2023growth_rates_1_LacI_TetRTup1_t_10000_aTc.mat')
load('27-May-2023growth_rates_1_LacI_TetRTup1_t_10000_IPTG.mat')
load('27-May-2023growth_rates_1_LacI_TetRTup1_t_10000.mat')
load('27-May-2023growth_rates_2_LacI_TetRTup1_t_10000_aTc.mat')
load('27-May-2023growth_rates_2_LacI_TetRTup1_t_10000_IPTG.mat')
load('27-May-2023growth_rates_2_LacI_TetRTup1_t_10000.mat')
load('27-May-2023growth_rates_3_LacI_TetRTup1_t_10000_aTc.mat')
load('27-May-2023growth_rates_3_LacI_TetRTup1_t_10000_IPTG.mat')
load('27-May-2023growth_rates_3_LacI_TetRTup1_t_10000.mat')
load('27-May-2023SimLacIfreeValues_aTc_rep1_combined.mat')
load('27-May-2023SimLacIfreeValues_aTc_rep2_combined.mat')
load('27-May-2023SimLacIfreeValues_aTc_rep3_combined.mat')
load('27-May-2023SimLacIfreeValues_IPTG_rep1_combined.mat')
load('27-May-2023SimLacIfreeValues_IPTG_rep2_combined.mat')
load('27-May-2023SimLacIfreeValues_IPTG_rep3_combined.mat')
load('27-May-2023SimLacIValues_aTc_rep1_combined.mat')
load('27-May-2023SimLacIValues_aTc_rep2_combined.mat')
load('27-May-2023SimLacIValues_aTc_rep3_combined.mat')
load('27-May-2023SimLacIValues_IPTG_rep1_combined.mat')
load('27-May-2023SimLacIValues_IPTG_rep2_combined.mat')
load('27-May-2023SimLacIValues_IPTG_rep3_combined.mat')
load('27-May-2023SimTetRTup1freeValues_IPTG_rep1_combined.mat')
load('27-May-2023SimTetRTup1freeValues_IPTG_rep2_combined.mat')
load('27-May-2023SimTetRTup1freeValues_IPTG_rep3_combined.mat')
load('27-May-2023SimTetRTup1Values_aTc_rep1.mat')
load('27-May-2023SimTetRTup1Values_aTc_rep2.mat')
load('27-May-2023SimTetRTup1Values_aTc_rep3.mat')
load('27-May-2023SimTetRTup1Values_IPTG_rep1_combined.mat')
load('27-May-2023SimTetRTup1Values_IPTG_rep2_combined.mat')
load('27-May-2023SimTetRTup1Values_IPTG_rep3_combined.mat')
load('27-May-2023time_aTc_rep1.mat')
load('27-May-2023time_aTc_rep2.mat')
load('27-May-2023time_aTc_rep3.mat')
load('27-May-2023time_IPTG_rep1.mat')
load('27-May-2023time_IPTG_rep2.mat')
load('27-May-2023time_IPTG_rep3.mat')

%Time dynamics for growth rate over aTc dose Rep 1
C = linspecer(height(data.dose));
for k = 1:width(SimFluoValues1_over_time_aTc_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep1_combined{a,k},SimFluoValues1_over_time_aTc_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 30)
        ylabel('mean Fluorescence','FontSize', 30)
        title('LacI Rep1 & TetRTup1 - aTc time dynamics', 'FontSize',40)
        
    
    end
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size            
    legend("show", 'Location', 'northeastoutside', 'FontSize',25)


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Time_dynamics/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep1_aTc_time_dynamics', '.jpg']);
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Time_dynamics/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep1_aTc_time_dynamics', '.fig']);

end
close all

%Time dynamics for growth rate over aTc dose Rep 2
C = linspecer(height(data.dose));
for k = 1:width(SimFluoValues2_over_time_aTc_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep2_combined{a,k},SimFluoValues2_over_time_aTc_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 30)
        ylabel('mean Fluorescence','FontSize', 30)
        title('LacI Rep2 & TetRTup1 - aTc time dynamics', 'FontSize',40)
        
    
    end
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size            
    legend("show", 'Location', 'northeastoutside', 'FontSize',25)


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Time_dynamics/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep1_aTc_time_dynamics', '.jpg']);
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Time_dynamics/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep1_aTc_time_dynamics', '.fig']);

end
close all

%Time dynamics for growth rate over IPTG dose Rep 1
%Time dynamics for growth rate over aTc dose Rep 2
C = linspecer(height(data.dose));
for k = 1:width(SimFluoValues3_over_time_aTc_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep3_combined{a,k},SimFluoValues3_over_time_aTc_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 30)
        ylabel('mean Fluorescence','FontSize', 30)
        title('LacI Rep3 & TetRTup1 - aTc time dynamics', 'FontSize',40)
        
    
    end
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size            
    legend("show", 'Location', 'northeastoutside', 'FontSize',25)


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Time_dynamics/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep1_aTc_time_dynamics', '.jpg']);
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Time_dynamics/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep1_aTc_time_dynamics', '.fig']);

end
close all

%Time dynamics for growth rate over IPTG dose Rep 1
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimFluoValues1_over_time_IPTG_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep1_combined{a,k},SimFluoValues1_over_time_IPTG_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 50)
        ylabel('Growth rate ','FontSize', 50)
        title('LacI Rep1 & TetRTup1 - IPTG time dynamics', 'FontSize',60)

    
    end
      ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size          
    %     legend("show", 'Location', 'northeastoutside', 'FontSize',25)
    colors = C; % yellow
    colormap(colors);

% Create colorbar
    cb = colorbar;

% Set colorbar ticks and labels
    cb.Ticks = [0 0.5 1]; % positions of the midpoints of the color segments
    cb.TickLabels = {data_IPTG.dose(1), data_IPTG.dose(50), data_IPTG.dose(100)};
    cb.FontSize = 25;
    ylabel(cb, 'IPTG dose','FontSize', 30)


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep2_IPTG_time_dynamics', '.jpg']);

end
close all

%Time dynamics for growth rate over IPTG dose Rep 3
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimFluoValues2_over_time_IPTG_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep2_combined{a,k},SimFluoValues2_over_time_IPTG_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 50)
        ylabel('Growth rate ','FontSize', 50)
        title('LacI Rep2 & TetRTup1 - IPTG time dynamics', 'FontSize',60)

    
    end
      ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size          
    %     legend("show", 'Location', 'northeastoutside', 'FontSize',25)
    colors = C; % yellow
    colormap(colors);

% Create colorbar
    cb = colorbar;

% Set colorbar ticks and labels
    cb.Ticks = [0 0.5 1]; % positions of the midpoints of the color segments
    cb.TickLabels = {data_IPTG.dose(1), data_IPTG.dose(50), data_IPTG.dose(100)};
    cb.FontSize = 25;
    ylabel(cb, 'IPTG dose','FontSize', 30)


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep2_IPTG_time_dynamics', '.jpg']);

end
close all


%Time dynamics for growth rate over IPTG dose Rep 3
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimFluoValues3_over_time_IPTG_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep3_combined{a,k},SimFluoValues3_over_time_IPTG_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 50)
        ylabel('Growth rate ','FontSize', 50)
        title('LacI Rep3 & TetRTup1 - IPTG time dynamics', 'FontSize',60)

    
    end
      ax = gca; % Get the current axes
    ax.XAxis.FontSize = 25; % X axis font size
    ax.YAxis.FontSize = 25; % Y axis font size          
    %     legend("show", 'Location', 'northeastoutside', 'FontSize',25)
    colors = C; % yellow
    colormap(colors);

% Create colorbar
    cb = colorbar;

% Set colorbar ticks and labels
    cb.Ticks = [0 0.5 1]; % positions of the midpoints of the color segments
    cb.TickLabels = {data_IPTG.dose(1), data_IPTG.dose(50), data_IPTG.dose(100)};
    cb.FontSize = 25;
    ylabel(cb, 'IPTG dose','FontSize', 30)


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_PMA1_rep2_IPTG_time_dynamics', '.jpg']);

end
close all

%Time dynamics for LacI over aTc dose Rep 1
C = linspecer(height(data.dose));
for k = 1:width(SimLacIValues_aTc_rep1_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep1_combined{a,k},SimLacIValues_aTc_rep1_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('LacI','FontSize', 18)
        title('LacI (1st repression coefficient) & TetRTup1 - aTc LacI time dynamics', 'FontSize',20)
        
    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_LacI_rep1_aTc_time_dynamics', '.jpg']);

end
close all

%Time dynamics for LacI over aTc dose Rep 2
C = linspecer(height(data.dose));
for k = 1:width(SimLacIValues_aTc_rep2_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep2_combined{a,k},SimLacIValues_aTc_rep2_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('LacI','FontSize', 18)
        title('LacI (2nd repression coefficient) & TetRTup1 - aTc LacI time dynamics', 'FontSize',20)
        
    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_LacI_rep2_aTc_time_dynamics', '.jpg']);

end
close all

%Time dynamics for LacI over aTc dose Rep 3
C = linspecer(height(data.dose));
for k = 1:width(SimLacIValues_aTc_rep3_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep3_combined{a,k},SimLacIValues_aTc_rep3_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('LacI','FontSize', 18)
        title('LacI (3rd repression coefficient) & TetRTup1 - aTc LacI time dynamics', 'FontSize',20)
        
    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_LacI_rep3_aTc_time_dynamics', '.jpg']);

end
close all

%Time dynamics for LacI over IPTG dose Rep 1
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimLacIValues_IPTG_rep1_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep1_combined{a,k},SimLacIValues_IPTG_rep1_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('LacI','FontSize', 18)
        title('LacI (1st repression coefficient) & TetRTup1 - IPTG LacI dynamics', 'FontSize',20)

    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_LacI_rep1_IPTG_time_dynamics', '.jpg']);

end
close all

%Time dynamics for LacI over IPTG dose Rep 2
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimLacIValues_IPTG_rep2_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep2_combined{a,k},SimLacIValues_IPTG_rep2_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('LacI','FontSize', 18)
        title('LacI (2nd repression coefficient) & TetRTup1 - IPTG LacI dynamics', 'FontSize',20)

    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_LacI_rep2_IPTG_time_dynamics', '.jpg']);

end
close all

%Time dynamics for LacI over IPTG dose Rep 3
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimLacIValues_IPTG_rep3_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep3_combined{a,k},SimLacIValues_IPTG_rep3_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('LacI','FontSize', 18)
        title('LacI (3rd repression coefficient) & TetRTup1 - IPTG LacI dynamics', 'FontSize',20)

    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_LacI_rep3_IPTG_time_dynamics', '.jpg']);

end
close all

 

%Time dynamics for LacIfree over IPTG dose Rep 1
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimLacIfreeValues_IPTG_rep1_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep1_combined{a,k},SimLacIfreeValues_IPTG_rep1_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('LacIfree','FontSize', 18)
        title('LacI (1st repression coefficient) & TetRTup1 - IPTG LacIfree dynamics', 'FontSize',20)

    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_LacIfree_rep1_IPTG_time_dynamics', '.jpg']);

end
close all

%Time dynamics for TetRTup1 over aTc dose Rep 1
C = linspecer(height(data.dose));
for k = 1:width(SimTetRTup1Values_aTc_rep1_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep1_combined{a,k},SimTetRTup1Values_aTc_rep1_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('TetRTup1','FontSize', 18)
        title('LacI (1st repression coefficient) & TetRTup1 - aTc TetRTup1 time dynamics', 'FontSize',20)
        
    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_TetRTup1_rep1_aTc_time_dynamics', '.jpg']);

end
close all

%Time dynamics for TetRTup1 over aTc dose Rep 2
C = linspecer(height(data.dose));
for k = 1:width(SimTetRTup1Values_aTc_rep2_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep2_combined{a,k},SimTetRTup1Values_aTc_rep2_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('TetRTup1','FontSize', 18)
        title('LacI (2nd repression coefficient) & TetRTup1 - aTc TetRTup1 time dynamics', 'FontSize',20)
        
    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_TetRTup1_rep2_aTc_time_dynamics', '.jpg']);

end
close all

%Time dynamics for TetRTup1 over aTc dose Rep 3
C = linspecer(height(data.dose));
for k = 1:width(SimTetRTup1Values_aTc_rep3_combined)
    figure(k);
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        % Convert the rounded value to a string and concatenate with the label
        plot(time_aTc_rep3_combined{a,k},SimTetRTup1Values_aTc_rep3_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('TetRTup1','FontSize', 18)
        title('LacI (3rd repression coefficient) & TetRTup1 - aTc TetRTup1 time dynamics', 'FontSize',20)
        
    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_TetRTup1_rep3_aTc_time_dynamics', '.jpg']);

end
close all

%Time dynamics for TetRTup1 over IPTG dose Rep 1
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimTetRTup1Values_IPTG_rep1_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep1_combined{a,k},SimTetRTup1Values_IPTG_rep1_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('TetRTup1','FontSize', 18)
        title('LacI (1st repression coefficient) & TetRTup1 - IPTG TetRTup1 dynamics', 'FontSize',20)

    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_TetRTup1_rep1_IPTG_time_dynamics', '.jpg']);

end
close all

%Time dynamics for TetRTup1 over IPTG dose Rep 2
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimTetRTup1Values_IPTG_rep2_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep2_combined{a,k},SimTetRTup1Values_IPTG_rep2_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('TetRTup1','FontSize', 18)
        title('LacI (2nd repression coefficient) & TetRTup1 - IPTG TetRTup1 dynamics', 'FontSize',20)

    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_TetRTup1_rep2_IPTG_time_dynamics', '.jpg']);

end
close all

%Time dynamics for TetRTup1 over IPTG dose Rep 3
C = linspecer(height(data_IPTG.dose));
for k = 1:width(SimTetRTup1Values_IPTG_rep3_combined)
    figure(k);
    for a = 1:height(data_IPTG.dose)
        hold on;

        rounded_dose = round(data_IPTG.dose(a) * 100) / 100;

        plot(time_IPTG_rep3_combined{a,k},SimTetRTup1Values_IPTG_rep3_combined{a,k},'-', 'LineWidth', 2, 'DisplayName', strcat('IPTG dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
        
        xlabel('Time', 'FontSize', 18)
        ylabel('TetRTup1','FontSize', 18)
        title('LacI (3rd repression coefficient) & TetRTup1 - IPTG TetRTup1 dynamics', 'FontSize',20)

    
    end
                
    legend("show", 'Location', 'northeastoutside')


%     set(figure(k), 'Position', get(0, 'Screensize'));
%     %save figure
%     saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy-'),num2str(k),'LacI_TetRTup1_TetRTup1_rep3_IPTG_time_dynamics', '.jpg']);

end
% close all


