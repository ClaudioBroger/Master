

p1 = [];
p2 = [];
p3 = [];
C = linspecer(height(data.dose));
for draw = 1:30


                figure(draw)
                for a = 1:height(data.dose)
                    hold on;
                    
                    % Round the data.dose(a) value to 2 decimal places
                    rounded_dose = round(data.dose(a) * 100) / 100;
%                     if draw == 1
                        SimFluoValues1 = table2array(SimFluoValues1combined(SimFluoValues1combined.Draw == draw,:));
                        % Convert the rounded value to a string and concatenate with the label
                        p1 =plot(log10(data_IPTG.dose),SimFluoValues1(a,1:100),'-', 'LineWidth', 2, 'Color', 'red');
                        
                        SimFluoValues2 = table2array(SimFluoValues2combined(SimFluoValues2combined.Draw == draw,:));
                        % Convert the rounded value to a string and concatenate with the label
                        p2 =plot(log10(data_IPTG.dose),SimFluoValues2(a,1:100),'-', 'LineWidth', 2, 'Color', 'blue');
                        
                         SimFluoValues3 = table2array(SimFluoValues3combined(SimFluoValues3combined.Draw == draw,:));
                        % Convert the rounded value to a string and concatenate with the label
                        p3 =plot(log10(data_IPTG.dose),SimFluoValues3(a,1:100),'-', 'LineWidth', 2, 'Color', 'green');
%                    
                    legend([p1 p2 p3],{'Repression coefficient 1', 'Repression coefficient 2', 'Repression coefficient 3'},'Location', 'northeastoutside','FontSize',25);
                    xlabel('log IPTG (nM)', 'FontSize', 60)
                    ylabel('mean Fluorescence','FontSize', 60)
                    title('Uncertainty all Rep - LacI TetRTup1 with aTc', 'FontSize',70)
                    ax = gca; % Get the current axes
                    ax.XAxis.FontSize = 40; % X axis font size
                    ax.YAxis.FontSize = 40; % Y axis font size  
                end

                    
                
                end


