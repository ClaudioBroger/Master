% Define your color vectors here
C1 = [1, 0, 0];  % Red for Fluovalues1
C2 = [0, 1, 0];  % Green for Fluovalues2
C3 = [0, 0, 1];  % Blue for Fluovalues3

for k = 1:height(rand_parameter)
    figure(k)
    for a = 1:height(data.dose)
        hold on;
        
        % Round the data.dose(a) value to 2 decimal places
        rounded_dose = round(data.dose(a) * 100) / 100;
        
        Fluovalues1 = SimFluoValues1combined(SimFluoValues1combined.Draw == k,:);
        Fluovalues1 = table2array(Fluovalues1(:,1:100));
        Fluovalues2 = SimFluoValues2combined(SimFluoValues2combined.Draw == k,:);
        Fluovalues2 = table2array(Fluovalues2(:,1:100));
        Fluovalues3 = SimFluoValues3combined(SimFluoValues3combined.Draw == k,:);
        Fluovalues3 = table2array(Fluovalues3(:,1:100));
        %LacIrep = linspace(min(SimFluoValues3combined.rep), max(SimFluoValues1combined.rep), 100);
        LacIrep1 = linspace(min(SimFluoValues1combined.rep), max(SimFluoValues1combined.rep), 100);
        LacIrep2 = linspace(min(SimFluoValues2combined.rep), max(SimFluoValues2combined.rep), 100);
        LacIrep3 = linspace(min(SimFluoValues3combined.rep), max(SimFluoValues3combined.rep), 100);

        
        % Convert the rounded value to a string and concatenate with the label
        % Use DisplayName for the first element of each group only
        if a == 1
            plot(LacIrep1,Fluovalues1(a,:),'-', 'LineWidth', 2, 'Color', C1, 'DisplayName', 'repression coefficient 1');
            plot(LacIrep2,Fluovalues2(a,:),'-', 'LineWidth', 2, 'Color', C2, 'DisplayName', 'repression coefficient 2');
            plot(LacIrep3,Fluovalues3(a,:),'-', 'LineWidth', 2, 'DisplayName', 'repression coefficient 3', 'Color', C3);
        else
            plot(LacIrep1,Fluovalues1(a,:),'-', 'LineWidth', 2, 'Color', C1, 'HandleVisibility', 'off');
            plot(LacIrep2,Fluovalues2(a,:),'-', 'LineWidth', 2, 'Color', C2, 'HandleVisibility', 'off');
            plot(LacIrep3,Fluovalues3(a,:),'-', 'LineWidth', 2, 'HandleVisibility', 'off', 'Color', C3);
        end
       
        xlabel('repression coefficient', 'FontSize', 18)
        ylabel('mean Fluorescence','FontSize', 18)
        %xlim([0.95 1.05]);
        title('LacI (all repression coefficients) & TetRTup1', 'FontSize',20)
        
    end

    legend("show", 'Location', 'northeastoutside')
end


            figure(k)
            for a = 1:height(data.dose)
                %subplot(3,1,2)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                Fluovalues2 = SimFluoValues2combined(SimFluoValues2combined.Draw == k,:);
                Fluovalues2 = table2array(Fluovalues2(:,1:100));
                %LacIrep = linspace(min(SimFluoValuescombined.rep), max(SimFluoValues3combined.rep), 100);

                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(LacIrep),Fluovalues2(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence','FontSize', 18)
                title('LacI (1st repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end

            


            figure(k)
            for a = 1:height(data.dose)
                %subplot(3,1,3)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                Fluovalues3 = SimFluoValues3combined(SimFluoValues3combined.Draw == k,:);
                Fluovalues3 = table2array(Fluovalues3(:,1:100));
                %LacIrep = linspace(min(SimFluoValues1combined.rep), max(SimFluoValues3combined.rep), 100);

                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(LacIrep),Fluovalues3(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence','FontSize', 18)
                title('LacI (1st repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end

            
            legend("show", 'Location', 'northeastoutside')
end