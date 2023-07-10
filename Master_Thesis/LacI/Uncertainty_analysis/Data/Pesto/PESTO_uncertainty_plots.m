%PESTO uncertainty

load('0324_PestoResults_cleanup.mat')

dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
%load model for uncertainty analysis
model = sbioloadproject('LacImodel_uncertainty');
FlPlot = true;
num_draws = 30;
counter = 0;

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));

parameters.name(1) = {"kLacI"};
parameters.name(5) = {"kCit"};
parameters.name(7) = {"dCit"};
parameters.name(8) = {"CitL"};
parameters.name(10) = {"LacIrep"};
parameters.name(11) = {"LacIrep2"};
parameters.name(12) = {"LacIrep3"};


parameters.name = cell2table(parameters.name);

[~,index] = ismember(ParaNames,parameters.name{:,:});
index = nonzeros(index);

min_par = parameters.min(index);
max_par = parameters.max(index);

CI = parameters.CI.S(index,:);

prob = parameters.S.par(index,:)';

for draw=1:num_draws
    for k=1:length(CI)
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
ParaNames(6:7) = [];
% for k = 1:num_draws
%     for j = 1:length(CI)
%         rand_parameter(k,j) = CI(j,1) + (CI(j,2) - CI(j,1)).*rand;
%     end
% end

rand_parameter(:,1:6) = 10.^rand_parameter(:,1:6);
rand_parameter(:,8:11) = 10.^rand_parameter(:,8:11);

rand_parameter = array2table(rand_parameter);
rand_parameter.Properties.VariableNames = ParaNames;


rand_parameter.mu(:) = 0.0077;
rand_parameter.kmaturation(:) = 0.0173;


rand_parameter.nMperUnit(:) = 10^1.1;
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
rand_parameter = rand_parameter(:,[string(ParaNames)]);
rand_parameter.dCit = zeros(num_draws,1);

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
            dataPos = strcat(dataPath, "data.mat");
            data1 = load(dataPos);
            data1.dose = [linspace(min(data1.dose), 10^3, 250) linspace(1+10^3, 10^7, 10000) linspace(1+10^7, max(data1.dose), 100)];
            data1.dose = data1.dose';
            
            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues1 = simulate_DR_IPTG(para,data1,ParaNames,model);
            DataMeans = data1.data.means;
            DataStd = data1.data.std;
            data1.dose(1) = 1;
            

%         %plot
        if FlPlot
            figure(2)
            hold on;
            plot(log10(data1.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 1st repression coefficient - parameter estimation with PESTO', 'FontSize',20)
            
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
            data2.dose = [linspace(min(data1.dose), 10^5, 100) linspace(1+10^5, 10^7, 1500) linspace(1+10^7, max(data1.dose), 100)];
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
            title('P4Lacn.2-cit + PAct1-LacI, 2nd repression coefficient - parameter estimation with PESTO', 'FontSize',20)
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
            data3.dose = [linspace(min(data1.dose), 10^4, 100) linspace(1+10^4, 10^9, 10000) linspace(1+10^9, 10^10, 100)];
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
            title('P4Lacn.2-cit + PAct1-LacI, 3rd repression coefficient - parameter estimation with PESTO', 'FontSize',20)
        end
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
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50_Pesto_30', '.mat'], 'SimFluoValues1combined');
%             writetable(SimFluoValues1combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_ec50', '.csv']);
%              %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50_Pesto_30', '.mat'], 'SimFluoValues2combined');
%             writetable(SimFluoValues2combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_ec50', '.csv']);
%             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_ec50_Pesto_30', '.mat'], 'SimFluoValues3combined');
%             writetable(SimFluoValues3combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_ec50', '.csv']);
%            
  


counter = 0;
for k = 1:num_draws
    x = SimFluoValues1combined.dose(SimFluoValues1combined.Draw == k);
    y = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_1_pesto(k,1) = xq(:,closest_index);
    ec50_1_pesto(k,2) = log10(ec50_1_pesto(k,1));
    ec50_1_pesto(k,3) = mid_fluo;
    ec50_1_pesto(k,4) = interpol(closest_index);
    ec50_1_pesto(k,5) = (abs(ec50_1_pesto(k,3) - ec50_1_pesto(k,4))) / mid_fluo;
    ec50_1_pesto(k,6) = k;
    counter = counter +1
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'ec50_1_random_pesto_200', '.mat'], 'ec50_1_pesto');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues2combined.dose(SimFluoValues2combined.Draw == k);
    y = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 1 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_2_pesto(k,1) = xq(:,closest_index);
    ec50_2_pesto(k,2) = log10(ec50_2_pesto(k,1));
    ec50_2_pesto(k,3) = mid_fluo;
    ec50_2_pesto(k,4) = interpol(closest_index);
    ec50_2_pesto(k,5) = (abs(ec50_2_pesto(k,3) - ec50_2_pesto(k,4))) / mid_fluo;
    ec50_2_pesto(k,6) = k;
    counter = counter +1
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'ec50_2_random_pesto_200', '.mat'], 'ec50_2_pesto');

counter = 0;
for k = 1:num_draws
    x = SimFluoValues3combined.dose(SimFluoValues3combined.Draw == k);
    y = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k);
    mid_fluo = abs(min(y) + max(y) * 0.5);
    xq = (min(x) : 7 : max(x));
    interpol = interp1(x,y,xq);
    [minvalue closest_index] = min(abs(interpol - mid_fluo));
    ec50_3_pesto(k,1) = xq(:,closest_index);
    ec50_3_pesto(k,2) = log10(ec50_3_pesto(k,1));
    ec50_3_pesto(k,3) = mid_fluo;
    ec50_3_pesto(k,4) = interpol(closest_index);
    ec50_3_pesto(k,5) = (abs(ec50_3_pesto(k,3) - ec50_3_pesto(k,4))) / mid_fluo;
    ec50_3_pesto(k,6) = k;
    counter = counter +1
end
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'ec50_3_random_pesto_200', '.mat'], 'ec50_3_pesto');

C = linspecer(num_draws);
figure(2)
hold on
title('P4Lacn.2-cit + PAct1-LacI, 1st repression coefficient - from PESTO', 'FontSize',20)
for k = 1:num_draws
    plot(ec50_1_pesto(k,2), ec50_1_pesto(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end


figure(3)
hold on
title('P4Lacn.2-cit + PAct1-LacI, 2nd repression coefficient - from PESTO', 'FontSize',20)
for k = 1:num_draws
    plot(ec50_2_pesto(k,2), ec50_2_pesto(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end

figure(4)
hold on
title('P4Lacn.2-cit + PAct1-LacI, 3rd repression coefficient - from PESTO', 'FontSize',20)
for k = 1:num_draws
    plot(ec50_3_pesto(k,2), ec50_3_pesto(k,4), '-k*','Color', C(k,:),'DisplayName',strcat('EC50-', num2str(k)),'MarkerSize',12)
end

set(figure(2), 'Position', get(0, 'Screensize'));
% % %             ax = gca(figure(1));
% % %             ax.ColorOrder = cmp;
% % %             %save figure
            saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_rep1_Pesto_ec50', '.jpg']);
% %            
% % %             %set the figure to be sceensize
            set(figure(3), 'Position', get(0, 'Screensize'));
% % %             ax = gca(figure(2));
% % %             ax.ColorOrder = cmp;
% % %             %save figure
            saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_rep2_Pesto_ec50', '.jpg']);
% %             
% % %             %set the figure to be sceensize
            set(figure(4), 'Position', get(0, 'Screensize'));
% % %             ax = gca(figure(3));
% % %             ax.ColorOrder = cmp;
% % %             %save figure
            saveas(figure(4), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_rep3_Pesto_ec50', '.jpg']);
% %     


C = linspecer(30);
for draw = 1:30
    
    SimFluoValues1 = SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == draw);
    figure(2)
    hold on;
    plot(log10(data1.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)), 'Color', C(draw,:));
    %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
    xlabel('log IPTG (nM)', 'FontSize', 18)
    ylabel('mean Fluorescence','FontSize', 18)
    title('P4Lacn.2-cit + PAct1-LacI, 1st repression coefficient - parameter estimation with PESTO', 'FontSize',20)
    

    hold off
    legend("show", 'Location', 'northeastoutside')

    SimFluoValues2 = SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == draw);
    

    figure(3)
    hold on;
    plot(log10(data2.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)), 'Color', C(draw,:));
    %errorbar(log10(data.dose),DataMeans,DataStd,'o','HandleVisibility','off');
    xlabel('log IPTG (nM)', 'FontSize', 18)
    ylabel('mean Fluorescence','FontSize', 18)
    title('P4Lacn.2-cit + PAct1-LacI, 2nd repression coefficient - parameter estimation with PESTO', 'FontSize',20)

    hold off
    legend("show", 'Location', 'northeastoutside')
    
    SimFluoValues3 = SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == draw);
    
    figure(4)
    hold on;
    plot(log10(data3.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)), 'Color', C(draw,:));
    %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
    xlim([0 10]);
    xlabel('log IPTG (nM)', 'FontSize', 18)
    ylabel('mean Fluorescence','FontSize', 18)
    title('P4Lacn.2-cit + PAct1-LacI, 3rd repression coefficient - parameter estimation with PESTO', 'FontSize',20)

    hold off
    legend("show", 'Location', 'northeastoutside')

end