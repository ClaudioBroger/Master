C = linspecer(num_draws);
counter = 0;

time_aTc_rep1_combined = cell(height(data.dose), height(rand_parameter));
time_aTc_rep2_combined = cell(height(data.dose), height(rand_parameter));
time_aTc_rep3_combined = cell(height(data.dose), height(rand_parameter));

time_IPTG_rep1_combined = cell(height(data_IPTG.dose), height(rand_parameter));
time_IPTG_rep2_combined = cell(height(data_IPTG.dose), height(rand_parameter));
time_IPTG_rep3_combined = cell(height(data_IPTG.dose), height(rand_parameter));

SimFluoValues1_over_time_aTc_combined = cell(height(data.dose), height(rand_parameter));
SimFluoValues2_over_time_aTc_combined = cell(height(data.dose), height(rand_parameter));
SimFluoValues3_over_time_aTc_combined = cell(height(data.dose), height(rand_parameter));

SimFluoValues1_over_time_IPTG_combined = cell(height(data_IPTG.dose), height(rand_parameter));
SimFluoValues2_over_time_IPTG_combined = cell(height(data_IPTG.dose), height(rand_parameter));
SimFluoValues3_over_time_IPTG_combined = cell(height(data_IPTG.dose), height(rand_parameter));

for draw = 1:num_draws
    
    %define parameter values as the values randomly drawn
     paraValues = table2array(rand_parameter(draw,:));
     paraValues = paraValues';
    %TFlMode = FlMode;
    
     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;
    
        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

     

            [SimFluoValues1, time_IPTG1, time_aTc1, SimFluoValues_over_time_aTc1, SimFluoValues_over_time_IPTG1] = simulate_DR_IPTG_TetR_aTc_more_output(para,data_IPTG,data,ParaNames,model);
           
            

%         %plot

            figure(2)
            hold on;
            plot(log10(data_IPTG.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('LacI (1st repression coefficient) & TetRTup1 - with aTc', 'FontSize',20)
            

        hold off
        legend("show", 'Location', 'northeastoutside')
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

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep2', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            [SimFluoValues2, time_IPTG2, time_aTc2, SimFluoValues_over_time_aTc2, SimFluoValues_over_time_IPTG2] = simulate_DR_IPTG_TetR_aTc_more_output(para,data_IPTG,data,ParaNames,model);

           
            
        %plot

            figure(3)
            hold on;
            plot(log10(data_IPTG.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o','HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('LacI (2nd repression coefficient) & TetRTup1 - with aTc', 'FontSize',20)

        hold off
        legend("show", 'Location', 'northeastoutside')

     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep3', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            [SimFluoValues3, time_IPTG3, time_aTc3, SimFluoValues_over_time_aTc3, SimFluoValues_over_time_IPTG3] = simulate_DR_IPTG_TetR_aTc_more_output(para,data_IPTG,data,ParaNames,model);

            
            

        %plot

            figure(4)
            hold on;
            plot(log10(data_IPTG.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('LacI (3rd repression coefficient) & TetRTup1 - with aTc', 'FontSize',20)

        hold off
        legend("show", 'Location', 'northeastoutside')

 

        %set SimFluoValues to table
        SimFluoValues1 = array2table(SimFluoValues1);
        SimFluoValues2 = array2table(SimFluoValues2);
        SimFluoValues3 = array2table(SimFluoValues3);
%         
%         %add a column to specify the repression coefficient
        SimFluoValues1.Dose = zeros(height(SimFluoValues1),1);
        SimFluoValues2.Dose = zeros(height(SimFluoValues2),1);
        SimFluoValues3.Dose = zeros(height(SimFluoValues3),1);
%         
         
        SimFluoValues1.draw = zeros(height(SimFluoValues1),1);
        SimFluoValues2.draw = zeros(height(SimFluoValues2),1);
        SimFluoValues3.draw = zeros(height(SimFluoValues3),1);
% 
%         %set column names
          
        SimFluoValues1.Dose = data_IPTG.dose;
        SimFluoValues2.Dose = data_IPTG.dose;
        SimFluoValues3.Dose = data_IPTG.dose;

        SimFluoValues1.Draw(:) = draw;
        SimFluoValues2.Draw(:) = draw;
        SimFluoValues3.Draw(:) = draw;

           
%         %combine all tables to one
        if draw == 1
            SimFluoValues1combined = SimFluoValues1;
            SimFluoValues2combined = SimFluoValues2;
            SimFluoValues3combined = SimFluoValues3;
            
        else
            SimFluoValues1combined = [SimFluoValues1combined; SimFluoValues1];
            SimFluoValues2combined = [SimFluoValues2combined; SimFluoValues2];
            SimFluoValues3combined = [SimFluoValues3combined; SimFluoValues3];
        end

        time_aTc_rep1_combined(:,draw) = time_aTc1;
        time_aTc_rep2_combined(:,draw) = time_aTc2;
        time_aTc_rep3_combined(:,draw) = time_aTc3;

        time_IPTG_rep1_combined(:,draw) = time_IPTG1;
        time_IPTG_rep2_combined(:,draw) = time_IPTG2;
        time_IPTG_rep3_combined(:,draw) = time_IPTG3;

        SimFluoValues1_over_time_aTc_combined(:,draw) = SimFluoValues_over_time_aTc1;
        SimFluoValues2_over_time_aTc_combined(:,draw) = SimFluoValues_over_time_aTc2;
        SimFluoValues3_over_time_aTc_combined(:,draw) = SimFluoValues_over_time_aTc3;
        
        SimFluoValues1_over_time_IPTG_combined(:,draw) = SimFluoValues_over_time_IPTG1;
        SimFluoValues2_over_time_IPTG_combined(:,draw) = SimFluoValues_over_time_IPTG2;
        SimFluoValues3_over_time_IPTG_combined(:,draw) = SimFluoValues_over_time_IPTG3;
        
% 
%         
%         
% 
        

end

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues1_LacI_TetRTup1_with_aTc', '.mat'], 'SimFluoValues1combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues2_LacI_TetRTup1_with_aTc', '.mat'], 'SimFluoValues2combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues3_LacI_TetRTup1_with_aTc', '.mat'], 'SimFluoValues3combined');

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues1_LacI_TetRTup1_over_time_aTc', '.mat'], 'SimFluoValues1_over_time_aTc_combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues2_LacI_TetRTup1_over_time_aTc', '.mat'], 'SimFluoValues2_over_time_aTc_combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues3_LacI_TetRTup1_over_time_aTc', '.mat'], 'SimFluoValues3_over_time_aTc_combined');

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues1_LacI_TetRTup1_over_time_IPTG', '.mat'], 'SimFluoValues1_over_time_IPTG_combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues2_LacI_TetRTup1_over_time_IPTG', '.mat'], 'SimFluoValues2_over_time_IPTG_combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues3_LacI_TetRTup1_over_time_IPTG', '.mat'], 'SimFluoValues3_over_time_IPTG_combined');
