time_IPTG1_combined = cell(height(data.dose), height(rand_parameter));
time_IPTG2_combined = cell(height(data.dose), height(rand_parameter));
time_IPTG3_combined = cell(height(data.dose), height(rand_parameter));

SimFluoValues_time1_combined = cell(height(data.dose), height(rand_parameter));
SimFluoValues_time2_combined = cell(height(data.dose), height(rand_parameter));
SimFluoValues_time3_combined = cell(height(data.dose), height(rand_parameter));

C = linspecer(num_draws);
counter = 0;
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

            [SimFluoValues1,  SimFluoValues_time1, time_IPTG1] = simulate_DR_IPTG_TetR_no_aTc(para, data,ParaNames,model);
           
            

%         %plot

            figure(2)
            hold on;
            plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)), 'Color', C(draw,:));

            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('LacI (1st repression coefficient) & TetRTup1', 'FontSize',20)
            

        hold off
        legend("show", 'Location', 'northeastoutside')

        

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

            [SimFluoValues2,  SimFluoValues_time2, time_IPTG2] = simulate_DR_IPTG_TetR_no_aTc(para,data,ParaNames,model);

           
            
        %plot

            figure(3)
            hold on;
            plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)), 'Color', C(draw,:));
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('LacI (2nd repression coefficient) & TetRTup1', 'FontSize',20)

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

            [SimFluoValues3,  SimFluoValues_time3, time_IPTG3] = simulate_DR_IPTG_TetR_no_aTc(para, data,ParaNames,model);



            

        %plot

            figure(4)
            hold on;
            plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)), 'Color', C(draw,:));
            xlim([0 10]);
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('LacI (3rd repression coefficient) & TetRTup1', 'FontSize',20)

        hold off
        legend("show", 'Location', 'northeastoutside')

 
%         
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
        SimFluoValues1.Repression_Coefficient = zeros(height(SimFluoValues1),1);
        SimFluoValues2.Repression_Coefficient = zeros(height(SimFluoValues2),1);
        SimFluoValues3.Repression_Coefficient = zeros(height(SimFluoValues3),1);
        
        SimFluoValues1.draw = zeros(height(SimFluoValues1),1);
        SimFluoValues2.draw = zeros(height(SimFluoValues2),1);
        SimFluoValues3.draw = zeros(height(SimFluoValues3),1);
% 
%         %set column names
        SimFluoValues1.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        SimFluoValues2.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        SimFluoValues3.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        
        SimFluoValues1.Dose = data.dose;
        SimFluoValues2.Dose = data.dose;
        SimFluoValues3.Dose = data.dose;

        SimFluoValues1.Draw(:) = draw;
        SimFluoValues2.Draw(:) = draw;
        SimFluoValues3.Draw(:) = draw;
% 
%         %convert second column to string 
        SimFluoValues1 = convertvars(SimFluoValues1, "Repression_Coefficient", 'string');
        SimFluoValues2 = convertvars(SimFluoValues2, "Repression_Coefficient", 'string');
        SimFluoValues3 = convertvars(SimFluoValues3, "Repression_Coefficient", 'string');
        
        %insert the appropriate repression coefficient
        SimFluoValues1.Repression_Coefficient(:) = 'LacIrep1';
        SimFluoValues2.Repression_Coefficient(:) = 'LacIrep2';
        SimFluoValues3.Repression_Coefficient(:) = 'LacIrep3';
%         
%         %combine all tables to one
        if draw == 1
            SimFluoValues1combined = SimFluoValues1;
            SimFluoValues2combined = SimFluoValues2;
            SimFluoValues3combined = SimFluoValues3;
            SimFluoValuescombined = [SimFluoValues1; SimFluoValues2; SimFluoValues3];
        else
            SimFluoValues1combined = [SimFluoValues1combined; SimFluoValues1];
            SimFluoValues2combined = [SimFluoValues2combined; SimFluoValues2];
            SimFluoValues3combined = [SimFluoValues3combined; SimFluoValues3];
            SimFluoValuescombined = [SimFluoValuescombined; SimFluoValues1; SimFluoValues2; SimFluoValues3];
        end


        time_IPTG1_combined(:,draw) = time_IPTG1;
        time_IPTG2_combined(:,draw) = time_IPTG2;
        time_IPTG3_combined(:,draw) = time_IPTG3;

        SimFluoValues_time1_combined(:,draw) = SimFluoValues_time1;
        SimFluoValues_time2_combined(:,draw) = SimFluoValues_time2;
        SimFluoValues_time3_combined(:,draw) = SimFluoValues_time3;


            counter = counter + 1


end

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues1_LacI_TetRTup1', '.mat'], 'SimFluoValues1combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues1_LacI_TetRTup2', '.mat'], 'SimFluoValues2combined');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/Parameter Mapping/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues1_LacI_TetRTup3', '.mat'], 'SimFluoValues3combined');