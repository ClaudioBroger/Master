C = linspecer(height(data.dose));
counter = 0;
% rand_parameter.kTetR = [];
% rand_parameter.dTetR = [];
for k = 1:height(rand_parameter)
    paraValues = table2array(rand_parameter(k,:));
     paraValues = paraValues';
    %TFlMode = FlMode;
    
     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;
    
        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep', 'p1', 'p2', 'growthMIN', 'growthMAX', 'f', 'g', 'kPMA1', 'dPMA1', 'PMAfactor', 'kLacTetRTup1'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

     

            growth_rate_rep1 = simulate_DR_IPTG_TetR_aTc_PMA1(para,data_IPTG,data,ParaNames,model);
           
%             if k == 1
%                 gr_rate_o_t_1 = growth_rate_over_time_rep1;
%             else
%                 gr_rate_o_t_1 = [gr_rate_o_t_1; growth_rate_over_time_rep1];
%             end
% 
%             figure(31)
%             subplot(3,1,1)
%             plot(time, growth_rate_over_time_rep1,'-', 'LineWidth', 2);
%                 
%             xlabel('Time', 'FontSize', 18)
%             ylabel('Growth rate','FontSize', 18)
%             title('LacI (1rd repression coefficient) & TetRTup1', 'FontSize',20)
            figure(k)
            for a = 1:height(data.dose)
                subplot(3,1,1)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(data_IPTG.dose),growth_rate_rep1(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('Growth rate','FontSize', 18)
                title('LacI (1st repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end
            
            legend("show", 'Location', 'northeastoutside')

%         %plot

            
            
        
        
%     
            % From here on the steps are repeated as in the module (with
            % first repression coefficient) above
        

        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep2', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit','mu', 'TetRTup1rep', 'p1', 'p2', 'growthMIN', 'growthMAX', 'f', 'g', 'kPMA1', 'dPMA1', 'PMAfactor', 'kLacTetRTup1'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            growth_rate_rep2 = simulate_DR_IPTG_TetR_aTc_PMA1(para,data_IPTG,data,ParaNames,model);
            
%             if k == 1
%                 gr_rate_o_t_2 = growth_rate_over_time_rep2;
%             else
%                 gr_rate_o_t_2 = [gr_rate_o_t_2; growth_rate_over_time_rep2];
%             end
% 
%             figure(31)
%             subplot(3,1,2)
%             plot(time, growth_rate_over_time_rep2,'-', 'LineWidth', 2);
%                 
%             xlabel('Time', 'FontSize', 18)
%             ylabel('Growth rate','FontSize', 18)
%             title('LacI (2rd repression coefficient) & TetRTup1', 'FontSize',20)
            
            figure(k)
            for a = 1:height(data.dose)
                subplot(3,1,2)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(data_IPTG.dose),growth_rate_rep2(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('Growth rate','FontSize', 18)
                title('LacI (2nd repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end
            
            legend("show", 'Location', 'northeastoutside')
           
            
        %plot

            

     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep3', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit','mu', 'TetRTup1rep', 'p1', 'p2', 'growthMIN', 'growthMAX', 'f', 'g', 'kPMA1', 'dPMA1', 'PMAfactor', 'kLacTetRTup1'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            growth_rate_rep3 = simulate_DR_IPTG_TetR_aTc_PMA1(para,data_IPTG,data,ParaNames,model);

%             if k == 1
%                 gr_rate_o_t_3 = growth_rate_over_time_rep3;
%             else
%                 gr_rate_o_t_3 = [gr_rate_o_t_3; growth_rate_over_time_rep3];
%             end
% 
%             figure(31)
%             subplot(3,1,3)
%             plot(time, growth_rate_over_time_rep3,'-', 'LineWidth', 2);
%                 
%             xlabel('Time', 'FontSize', 18)
%             ylabel('Growth rate','FontSize', 18)
%             title('LacI (3rd repression coefficient) & TetRTup1', 'FontSize',20)
        %plot

            figure(k)
            for a = 1:height(data.dose)
                subplot(3,1,3)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(data_IPTG.dose),growth_rate_rep3(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('Growth rate','FontSize', 18)
                title('LacI (3rd repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end
            
            legend("show", 'Location', 'northeastoutside')
            objf  = [];
            objfn = {};
            ndata = 0;

%         if FlPlot
%             figure(draw+3)
%             hold on;
%             plot(log10(data1.dose), SimFluoValues1, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)));
%             plot(log10(data2.dose), SimFluoValues2, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)));
%             plot(log10(data3.dose), SimFluoValues3, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)));
%             %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
%             xlabel('log IPTG (nM)', 'FontSize', 18)
%             ylabel('mean Fluorescence','FontSize', 18)
%             title('P4Lacn.2-cit + PAct1-LacI, all repression coefficients, draw',strcat(num2str(draw)), 'FontSize',20);
%         end
%         hold off
%         legend("show", 'Location', 'northeastoutside')
        %set SimFluoValues to table
        growth_rate_rep1 = array2table(growth_rate_rep1);
        growth_rate_rep2 = array2table(growth_rate_rep2);
        growth_rate_rep3 = array2table(growth_rate_rep3);
%         
%         %add a column to specify the repression coefficient
%         SimFluoValues1.Dose = zeros(height(SimFluoValues1),1);
%         SimFluoValues2.Dose = zeros(height(SimFluoValues2),1);
%         SimFluoValues3.Dose = zeros(height(SimFluoValues3),1);
% %         
%         SimFluoValues1.Repression_Coefficient = zeros(height(SimFluoValues1),1);
%         SimFluoValues2.Repression_Coefficient = zeros(height(SimFluoValues2),1);
%         SimFluoValues3.Repression_Coefficient = zeros(height(SimFluoValues3),1);
%         
%         SimFluoValues1.draw = zeros(height(SimFluoValues1),1);
%         SimFluoValues2.draw = zeros(height(SimFluoValues2),1);
%         SimFluoValues3.draw = zeros(height(SimFluoValues3),1);
% 
%         %set column names
         

        growth_rate_rep1 = growth_rate_rep1(1:height(data.dose),:);
        growth_rate_rep2 = growth_rate_rep2(1:height(data.dose),:);
        growth_rate_rep3 = growth_rate_rep3(1:height(data.dose),:);

        growth_rate_rep1.Dose = data.dose;
        growth_rate_rep2.Dose = data.dose;
        growth_rate_rep3.Dose = data.dose;

        growth_rate_rep1.Draw(:) = k;
        growth_rate_rep2.Draw(:) = k;
        growth_rate_rep3.Draw(:) = k;

        growth_rate_rep1.rep(:) = table2array(rand_parameter(k,"LacIrep"));
        growth_rate_rep2.rep(:) = table2array(rand_parameter(k,"LacIrep2"));
        growth_rate_rep3.rep(:) = table2array(rand_parameter(k,"LacIrep3"));
% 
%         %convert second column to string 
          
%         %combine all tables to one
        if k == 1
            growth_rate_rep1_combined = growth_rate_rep1;
            growth_rate_rep2_combined = growth_rate_rep2;
            growth_rate_rep3_combined = growth_rate_rep3;
           
        else
            growth_rate_rep1_combined = [growth_rate_rep1_combined; growth_rate_rep1];
            growth_rate_rep2_combined = [growth_rate_rep2_combined; growth_rate_rep2];
            growth_rate_rep3_combined = [growth_rate_rep3_combined; growth_rate_rep3];
            
        end
% 
%         
%         

        %set the figure to be sceensize
             set(figure(k), 'Position', get(0, 'Screensize'));
            %save figure
             saveas(figure(k), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Plots/Full_aTc/More_time/More_aTc/Corrected/', datestr(now, 'dd-mmm-yyyy'),num2str(k),'LacI_TetRTup1_PMA1_more_aTc_more_time', '.jpg']);

            counter = counter + 1


end

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Results/', datestr(now, 'dd-mmm-yyyy'),'growth_rates_1_LacI_TetRTup1_more_time_more_aTc', '.mat'], 'growth_rate_rep1_combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Results/', datestr(now, 'dd-mmm-yyyy'),'growth_rates_2_LacI_TetRTup1_more_time_more_aTc', '.mat'], 'growth_rate_rep2_combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/Results/', datestr(now, 'dd-mmm-yyyy'),'growth_rates_3_LacI_TetRTup1_more_time_more_aTc', '.mat'], 'growth_rate_rep3_combined');


