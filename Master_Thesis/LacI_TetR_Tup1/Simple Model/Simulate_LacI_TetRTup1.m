%PESTO uncertainty

load('1704_fitModel_Pesto.mat')

dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
dataPath_aTc = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/');
%load model for uncertainty analysis
model = sbioloadproject('LacI_TetR_Tup1_model');
FlPlot = true;
num_draws = 30;
counter = 0;

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));


parameters.name = cell2table(parameters.name);

[~,index] = ismember(ParaNames,parameters.name{:,:});
index = nonzeros(index);

min_par = parameters.min(index);
max_par = parameters.max(index);


prob = parameters.S.par(index,:)';

for draw=1:num_draws
    for k=1:length(min_par)
        pd = fitdist(prob(:,k),'Kernel');
         X = min_par(k):0.001:max_par(k);
         Y = pdf(pd,X);
%          figure(1)
%          subplot(2,6,k)
         plot(Y)
         rand_parameter(draw,k) = randsample(X,1,true,Y);
         
    end
end



% num_draws = 30;
% for k = 1:num_draws
%     for j = 1:length(CI)
%         rand_parameter(k,j) = CI(j,1) + (CI(j,2) - CI(j,1)).*rand;
%     end
% end



rand_parameter = array2table(rand_parameter);
rand_parameter.Properties.VariableNames = table2array(parameters.name);

rand_parameter.mu(:) = 0.0077;
rand_parameter.kmaturation(:) = 0.0173;


rand_parameter.nMperUnit(:) = 10^1.1;
rand_parameter.a(:) = 1;
rand_parameter.f(:) = 1;
rand_parameter.g(:) = 1;
rand_parameter.b(:) = 1;
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
rand_parameter = rand_parameter(:,[string(ParaNames)]);
rand_parameter.dCit(:) = 0;

rand_parameter(:,1:5) = array2table(10.^table2array(rand_parameter(:,1:5)));
rand_parameter(:,9:32) = array2table(10.^table2array(rand_parameter(:,9:32)));
rand_parameter.p1(:) = 5;
rand_parameter.p2(:) = 5;
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/',datestr(now, 'dd-mmm-yyyy'),'_rand_parameter_Pesto.mat'],'rand_parameter')

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
            dataPos_IPTG = strcat(dataPath, "data.mat");
            data_IPTG = load(dataPos_IPTG);
            data_IPTG.dose = [linspace(min(data_IPTG.dose), 10^3, 250) linspace(1+10^3, 10^7, 10000) linspace(1+10^7, max(data_IPTG.dose), 100)];
            data_IPTG.dose = data_IPTG.dose';
            dataPos_aTc = strcat(dataPath_aTc, "dataSRtup1.mat");
            data_aTc = load(dataPos_aTc);
            
            
            
            NamestoZero = setdiff(ParaNames,{'kactTetR','k7tetCit', 'kL7tet', 'dTetR', 'dCit', 'thetaTetR', 'nTetR', 'f', 'nTup1', 'KdTetR', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime' , 'thetaTup1', 'thetaLacI3mut_2', 'a', 'thetaLacI3mut', 'p1', 'p2'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues1 = simulate_DR_IPTG_aTc(para,data_IPTG,data_aTc,ParaNames,model);
            DataMeans = data_IPTG.data.means;
            DataStd = data_IPTG.data.std;
            data_IPTG.dose(1) = 1;
            data_aTc.dose(1) = 1;
            

%         %plot
        if FlPlot
            figure(2)
            hold on;
            plot(log10(data_IPTG.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 1st repression coefficient', 'FontSize',20)
            
        end
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
            dataPos = strcat(dataPath, "data.mat");
            data2 = load(dataPos);
            data2.dose = [linspace(min(data_IPTG.dose), 10^5, 100) linspace(1+10^5, 10^7, 1500) linspace(1+10^7, max(data_IPTG.dose), 100)];
            data2.dose = data2.dose';

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep2', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues2 = simulate_DR_IPTG(para,data2,ParaNames,model);
            DataMeans = data2.data.means;
            DataStd = data2.data.std;
            data2.dose(1) = 1;

           
            
        %plot
        if FlPlot
            figure(3)
            hold on;
            plot(log10(data2.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o','HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 2nd repression coefficient', 'FontSize',20)
        end
        hold off
        legend("show", 'Location', 'northeastoutside')

     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            dataPos = strcat(dataPath, "data.mat");
            data3 = load(dataPos);
            data3.dose = [linspace(min(data_IPTG.dose), 10^4, 100) linspace(1+10^4, 10^9, 10000) linspace(1+10^9, 10^10, 100)];
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
        if FlPlot
            figure(4)
            hold on;
            plot(log10(data3.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlim([0 10]);
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 3rd repression coefficient', 'FontSize',20)
        end
        hold off
        legend("show", 'Location', 'northeastoutside')

 
        if FlPlot
            figure(draw+3)
            hold on;
            plot(log10(data1.dose), SimFluoValues1, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)));
            plot(log10(data2.dose), SimFluoValues2, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)));
            plot(log10(data3.dose), SimFluoValues3, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, all repression coefficients, draw',strcat(num2str(draw)), 'FontSize',20);
        end
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
        
        SimFluoValues1.dose = data_IPTG.dose;
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
            set(figure(draw+3), 'Position', get(0, 'Screensize'));
            %save figure
            saveas(figure(draw+3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),num2str(draw),'Figure_draw_allrep', '.jpg']);

            counter = counter + 1
            fprintf('End run', counter)

end
%             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50_Pesto_200', '.mat'], 'SimFluoValues1combined');
%             writetable(SimFluoValues1combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50', '.csv']);
%              %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50_Pesto_200', '.mat'], 'SimFluoValues2combined');
%             writetable(SimFluoValues2combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50', '.csv']);
%             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_ec50_Pesto_200', '.mat'], 'SimFluoValues3combined');
%  