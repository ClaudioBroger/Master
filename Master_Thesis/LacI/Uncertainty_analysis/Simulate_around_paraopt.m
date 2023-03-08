%Load data
load('2023_02_09_fitModel_Hyperspace_npoints20_Eline.mat')
%Load model
dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
model = sbioloadproject('LacImodel');
FlPlot = true;
%Load model settings
ModelSettings_4;

[mincost,idx_mincost] = min(costlist);
paraopt = paraoptlist(idx_mincost,:);
new_names = ["kLacI", "kCit", "dLacI", "LacIrep", "KdLacI", "nLacI", "nMperUnit", "LacIrep2", "kLacI2", "CitL", "LacIrep3", "kLacI3", "kCit2", "kCit2L", "dLacI3", "nLacI2", "LacIrep3.1"];

%paraopt([1 2 3 4 5 6 8 9 10 11 12 13 14 15 17]) = log10(paraopt([1 2 3 4 5 6 8 9 10 11 12 13 14 15 17]));

min_para = 0.9*paraopt;
max_para = 1.1*paraopt;

%min_para([1 2 3 4 5 6 8 9 10 11 12 13 14 15 17]) = 10.^(min_para([1 2 3 4 5 6 8 9 10 11 12 13 14 15 17]));
%max_para([1 2 3 4 5 6 8 9 10 11 12 13 14 15 17]) = 10.^(max_para([1 2 3 4 5 6 8 9 10 11 12 13 14 15 17]));

num_draws = 20;

for k = 1:length(paraopt)
    for j = 1:num_draws
        range = min_para(k):0.000001:max_para(k);
        number = randi(length(min_para(k):0.000001:max_para(k)));
        paraopt_range(j,k) = range(number);
    end
end

paraopt_range = array2table(paraopt_range);
paraopt_range.Properties.VariableNames = new_names;

paraopt_range.dCit = zeros(num_draws,1);
for n = 1:num_draws
    paraopt_range.mu(n) = 0.0077;
    paraopt_range.kmaturation(n) = 0.0173;
end

model = sbioloadproject('LacImodel_uncertainty');

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
paraopt_range = paraopt_range(:,[string(ParaNames)]);

C = linspecer(num_draws);
for draw = 1:num_draws
    
     paraValues = table2array(paraopt_range(draw,:));
     paraValues = paraValues';
    %TFlMode = FlMode;
    
     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;
    %Simulate every module three times; each time with another repression
    %coefficient -> 15 simulations per random parameter draw
        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            dataPos = strcat(dataPath, "data.mat");
            data = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues1 = simulate_DR_IPTG(para,data,ParaNames,model);
            DataMeans = data.data.means;
            DataStd = data.data.std;
            data.dose(1) = 1;
            
        %plot
        if FlPlot
            figure(1)
            hold on;
            plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)),'Color',C(draw,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 1st repression coefficient', 'FontSize',20)
            
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
            data = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep2', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues2 = simulate_DR_IPTG(para,data,ParaNames,model);
            DataMeans = data.data.means;
            DataStd = data.data.std;
            data.dose(1) = 1;
            
        %plot
        if FlPlot
            figure(2)
            hold on;
            plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)),'Color',C(draw,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o','HandleVisibility','off');
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
            data = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues3 = simulate_DR_IPTG(para,data,ParaNames,model);
            DataMeans = data.data.means;
            DataStd = data.data.std;
            data.dose(1) = 1;
            
        %plot
        if FlPlot
            figure(3)
            hold on;
            plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)),'Color',C(draw,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 3rd repression coefficient', 'FontSize',20)
        end
        hold off
        legend("show", 'Location', 'northeastoutside')

        if FlPlot
            figure(draw+3)
            hold on;
            plot(log10(data.dose), SimFluoValues1, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)));
            plot(log10(data.dose), SimFluoValues2, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep2-', num2str(draw)));
            plot(log10(data.dose), SimFluoValues3, '-', 'LineWidth',2, 'DisplayName',strcat('Simulation-rep3-', num2str(draw)));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
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
        
        %add a column to specify the repression coefficient
        SimFluoValues1.Dose = zeros(height(SimFluoValues1),1);
        SimFluoValues2.Dose = zeros(height(SimFluoValues2),1);
        SimFluoValues3.Dose = zeros(height(SimFluoValues3),1);
        
        SimFluoValues1.Repression_Coefficient = zeros(height(SimFluoValues1),1);
        SimFluoValues2.Repression_Coefficient = zeros(height(SimFluoValues2),1);
        SimFluoValues3.Repression_Coefficient = zeros(height(SimFluoValues3),1);
        
        SimFluoValues1.draw = zeros(height(SimFluoValues1),1);
        SimFluoValues2.draw = zeros(height(SimFluoValues2),1);
        SimFluoValues3.draw = zeros(height(SimFluoValues3),1);

        %set column names
        SimFluoValues1.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        SimFluoValues2.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        SimFluoValues3.Properties.VariableNames = {'SimFluoValues', 'Dose', 'Repression_Coefficient', 'Draw'};
        
        SimFluoValues1.Dose = data.dose;
        SimFluoValues2.Dose = data.dose;
        SimFluoValues3.Dose = data.dose;

        SimFluoValues1.Draw(:) = draw;
        SimFluoValues2.Draw(:) = draw;
        SimFluoValues3.Draw(:) = draw;

        %convert second column to string 
        SimFluoValues1 = convertvars(SimFluoValues1, "Repression_Coefficient", 'string');
        SimFluoValues2 = convertvars(SimFluoValues2, "Repression_Coefficient", 'string');
        SimFluoValues3 = convertvars(SimFluoValues3, "Repression_Coefficient", 'string');
        
        %insert the appropriate repression coefficient
        SimFluoValues1.Repression_Coefficient(:) = 'LacIrep1';
        SimFluoValues2.Repression_Coefficient(:) = 'LacIrep2';
        SimFluoValues3.Repression_Coefficient(:) = 'LacIrep3';
        
        %combine all tables to one
        if draw == 1
            SimFluoValues1combined = SimFluoValues1;
            SimFluoValues2combined = SimFluoValues2;
            SimFluoValues3combined = SimFluoValues3;
        else
            SimFluoValues1combined = [SimFluoValues1combined; SimFluoValues1];
            SimFluoValues2combined = [SimFluoValues2combined; SimFluoValues2];
            SimFluoValues3combined = [SimFluoValues3combined; SimFluoValues3];
        end

        SimFluoValuescombined = [SimFluoValues1; SimFluoValues2; SimFluoValues3];
        

        %set the figure to be sceensize
            set(figure(draw+3), 'Position', get(0, 'Screensize'));
            %save figure
            saveas(figure(draw+3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),num2str(draw),'Figure_around_paraopt_draw_allrep', '.jpg']);
            %savefig(fig, '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_', num2str(draw))
            %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'SimulationFluoValues_around_paraopt_allrep', '.mat'], 'SimFluoValuescombined');


end
            %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_around_paraopt_rep1_', '.mat'], 'SimFluoValues1combined');
            writetable(SimFluoValues1combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep1_', '.csv']);
             %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_around_paraopt_rep2_', '.mat'], 'SimFluoValues2combined');
            writetable(SimFluoValues2combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep2_', '.csv']);
            %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_around_paraopt_rep3_', '.mat'], 'SimFluoValues3combined');
            writetable(SimFluoValues3combined,['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'SimulationFluoValues_rep3_', '.csv']);
           
            %func = @(x) colorspace('RGB->Lab',x);
            %cmp=distinguishable_colors(30);
            %set the figure to be sceensize
            set(figure(1), 'Position', get(0, 'Screensize'));
% %          ax = gca(figure(1));
% %          ax.ColorOrder = cmp;
            %save figure
            saveas(figure(1), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_around_paraopt_rep1_', '.jpg']);
           
            %set the figure to be sceensize
            set(figure(2), 'Position', get(0, 'Screensize'));
%             ax = gca(figure(2));
%             ax.ColorOrder = cmp;
            %save figure
            saveas(figure(2), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_around_paraopt_rep2_', '.jpg']);
            
            %set the figure to be sceensize
            set(figure(3), 'Position', get(0, 'Screensize'));
%             ax = gca(figure(3));
%             ax.ColorOrder = cmp;
            %save figure
            saveas(figure(3), ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/', datestr(now, 'dd-mmm-yyyy'),'_draw_', num2str(draw),'_Figure_around_paraopt_rep3_', '.jpg']);
           