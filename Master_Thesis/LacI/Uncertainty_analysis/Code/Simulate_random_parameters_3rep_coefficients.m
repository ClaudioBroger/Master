%set datapath to location of experimental data
dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
%load model for uncertainty analysis
model = sbioloadproject('LacImodel_uncertainty');
FlPlot = false;
num_draws = 200;
counter = 0;
% path = '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/';
% dataPos = strcat(path, "26-Mar-2023_rand_parameter_uncertainty_analysis_ec50.mat");
% load(dataPos)

%random parameters are in log scale -> need to scale them back to use them
%as parameters in the model
%hill coefficients are not in log scale -> not scale them back
rand_parameter = 10.^(table2array(rand_parameter));
rand_parameter(:,6) = log10(rand_parameter(:,6));
rand_parameter(:,16) = log10(rand_parameter(:,16));
rand_parameter = array2table(rand_parameter);
new_names = ["kLacI", "kCit", "dLacI", "LacIrep", "KdLacI", "nLacI", "nMperUnit", "LacIrep2", "kLacI2", "CitL", "LacIrep3", "kLacI3", "kCit2", "kCit2L", "dLacI3", "nLacI2", "LacIrep3.1"];
rand_parameter.Properties.VariableNames = new_names;
% 
% %add dCit as 0s because in the optimal parameter set is is set to 0
 rand_parameter.dCit = zeros(num_draws,1);
% %add constant values for mu and kmaturation
for n = 1:num_draws
    rand_parameter.mu(n) = 0.0077;
    rand_parameter.kmaturation(n) = 0.0173;
end

% SimFluoValues1combined = zeros(12, 4);
% SimFluoValues2combined = zeros(12, 4);
% SimFluoValues3combined = zeros(12, 4);
% 
% SimFluoValues1combined = array2table(SimFluoValues1combined);
% SimFluoValues2combined = array2table(SimFluoValues2combined);
% SimFluoValues3combined = array2table(SimFluoValues3combined);
% 
% SimFluoValues1combined.Properties.VariableNames = ["SimFluoValues", "Dose", "Repression_Coefficient", "Draw"];
% SimFluoValues2combined.Properties.VariableNames = ["SimFluoValues", "Dose", "Repression_Coefficient", "Draw"];
% SimFluoValues3combined.Properties.VariableNames = ["SimFluoValues", "Dose", "Repression_Coefficient", "Draw"];


%get the names of the parameters from the model for uncertainty analysis
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
rand_parameter = rand_parameter(:,[string(ParaNames)]);

%  save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/',datestr(now, 'dd-mmm-yyyy'),'_rand_parameter_uncertainty_analysis_ec50.mat'],'rand_parameter')

C = linspecer(num_draws);
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
            dataPos = strcat(dataPath, "data.mat");
            data1 = load(dataPos);
            data1.dose = [linspace(min(data1.dose), 10^3, 1000) linspace(1+10^3, 10^7, 10000) linspace(1+10^7, max(data1.dose), 1000)];
            data1.dose = data1.dose';
            
            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues1 = simulate_DR_IPTG(para,data1,ParaNames,model);
            DataMeans = data1.data.means;
            DataStd = data1.data.std;
            data1.dose(1) = 1;
            

        %plot
%         if FlPlot
%             figure(1)
%             hold on;
%             plot(log10(data1.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)), 'Color', C(draw,:));
%             %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
%             xlabel('log IPTG (nM)', 'FontSize', 18)
%             ylabel('mean Fluorescence','FontSize', 18)
%             title('P4Lacn.2-cit + PAct1-LacI, 1st repression coefficient', 'FontSize',20)
%             
%         end
%         hold off
%         legend("show", 'Location', 'northeastoutside')
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
            dataPos = strcat(dataPath, "data.mat");
            data2 = load(dataPos);
            data2.dose = [linspace(min(data2.dose), 10^5, 100) linspace(1+10^5, 10^7, 1500) linspace(1+10^7, max(data2.dose), 100)];
            data2.dose = data2.dose';

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep2', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues2 = simulate_DR_IPTG(para,data2,ParaNames,model);
            DataMeans = data2.data.means;
            DataStd = data2.data.std;
            data2.dose(1) = 1;

           
            
        %plot
%         if FlPlot
%             figure(2)
%             hold on;
%             plot(log10(data2.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)), 'Color', C(draw,:));
%             %errorbar(log10(data.dose),DataMeans,DataStd,'o','HandleVisibility','off');
%             xlabel('log IPTG (nM)', 'FontSize', 18)
%             ylabel('mean Fluorescence','FontSize', 18)
%             title('P4Lacn.2-cit + PAct1-LacI, 2nd repression coefficient', 'FontSize',20)
%         end
%         hold off
%         legend("show", 'Location', 'northeastoutside')

     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            dataPos = strcat(dataPath, "data.mat");
            data3 = load(dataPos);
            data3.dose = [linspace(min(data3.dose), 10^5.6, 100) linspace(1+10^6.5, 10^8.5, 1000) linspace(1+10^8.5, 10^10, 100)];
            data3.dose = data3.dose';

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues3 = simulate_DR_IPTG(para,data3,ParaNames,model);
            DataMeans = data3.data.means;
            DataStd = data3.data.std;
            data3.dose(1) = 1;

            %[hill3(draw,1) ec50_3(draw,1)] = doseResponse(data3.dose, SimFluoValues3);

            

        %plot
%         if FlPlot
%             figure(3)
%             hold on;
%             plot(log10(data3.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)), 'Color', C(draw,:));
%             %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
%             xlim([0 10]);
%             xlabel('log IPTG (nM)', 'FontSize', 18)
%             ylabel('mean Fluorescence','FontSize', 18)
%             title('P4Lacn.2-cit + PAct1-LacI, 3rd repression coefficient', 'FontSize',20)
%         end
%         hold off
%         legend("show", 'Location', 'northeastoutside')

 
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
        
        SimFluoValues1.dose = data1.dose;
        SimFluoValues2.dose = data2.dose;
        SimFluoValues3.dose = data3.dose;

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
% 
%         
%         
% 
%         %set the figure to be sceensize
%             set(figure(draw+3), 'Position', get(0, 'Screensize'));
%             %save figure
%             saveas(figure(draw+3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),num2str(draw),'Figure_draw_allrep', '.jpg']);

            counter = counter + 1
            fprintf('End run', counter)

end
%             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/30. March/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50_200', '.mat'], 'SimFluoValues1combined');
            writetable(SimFluoValues1combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/28. March/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50_200', '.csv']);
%              %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/30. March/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50_200', '.mat'], 'SimFluoValues2combined');
            writetable(SimFluoValues2combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/28. March/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50_200', '.csv']);
%             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/30. March/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_ec50_200', '.mat'], 'SimFluoValues3combined');
            writetable(SimFluoValues3combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/28. March/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_ec50_200', '.csv']);
%            
%             save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'SimulationFluoValues_allrep_ec50', '.mat'], 'SimFluoValuescombined');

            %func = @(x) colorspace('RGB->Lab',x);
            %cmp=distinguishable_colors(30);
            %set the figure to be sceensize
%             set(figure(1), 'Position', get(0, 'Screensize'));
% % %             ax = gca(figure(1));
% % %             ax.ColorOrder = cmp;
% % %             %save figure
%             saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_rep1_', '.jpg']);
%            
% % %             %set the figure to be sceensize
%             set(figure(2), 'Position', get(0, 'Screensize'));
% % %             ax = gca(figure(2));
% % %             ax.ColorOrder = cmp;
% % %             %save figure
%             saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_rep2_', '.jpg']);
%             
% % %             %set the figure to be sceensize
%             set(figure(3), 'Position', get(0, 'Screensize'));
% % %             ax = gca(figure(3));
% % %             ax.ColorOrder = cmp;
% % %             %save figure
%             saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_rep3_', '.jpg']);
%            
% 


%  

