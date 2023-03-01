%set datapath to location of experimental data
dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
model = sbioloadproject('LacImodel_uncertainty');
FlPlot = true;

rand_parameter = 10.^(table2array(rand_parameter));
rand_parameter(:,6) = log10(rand_parameter(:,6));
rand_parameter(:,16) = log10(rand_parameter(:,16));
rand_parameter = array2table(rand_parameter);
new_names = ["kLacI", "kCit", "dLacI", "LacIrep", "KdLacI", "nLacI", "nMperUnit", "LacIrep2", "kLacI2", "CitL", "LacIrep3", "kLacI3", "kCit2", "kCit2L", "dLacI3", "nLacI2", "LacIrep3.1"];
rand_parameter.Properties.VariableNames = new_names;

rand_parameter.dCit = zeros(num_draws,1);
for n = 1:num_draws
    rand_parameter.mu(n) = 0.0077;
    rand_parameter.kmaturation(n) = 0.0173;
end

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
rand_parameter = rand_parameter(:,[string(ParaNames)]);

for draw = 1:num_draws
    
     paraValues = table2array(rand_parameter(draw,:));
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

            SimFluoValues = simulate_DR_IPTG(para,data,ParaNames,model);
            DataMeans = data.data.means;
            DataStd = data.data.std;
            data.dose(1) = 1;
            
        %plot
        if FlPlot
            figure(1)
            subplot(4,3,draw);
            hold on;
            plot(log10(data.dose),SimFluoValues,'-', 'LineWidth', 2);
            errorbar(log10(data.dose),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI', 'FontSize',20)
            legend({'simulation'}, 'Location', 'southeast')
        end
    
    
            % From here on the steps are repeated as in the module (with
            % first repression coefficient) above
           
            %set the figure to be sceensize
            set(gcf, 'Position', get(0, 'Screensize'));
            %save figure
            saveas(gcf, ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_rep1_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.jpg']);
            %savefig(fig, '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_', num2str(draw))
            %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/SimulationFluoValues_draw_rep1_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.mat'], 'SimFluoValues');

end

for draw = 1:num_draws
    
     paraValues = table2array(rand_parameter(draw,:));
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

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep2', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues = simulate_DR_IPTG(para,data,ParaNames,model);
            DataMeans = data.data.means;
            DataStd = data.data.std;
            data.dose(1) = 1;
            
        %plot
        if FlPlot
            figure(2)
            subplot(4,3,draw);
            hold on;
            plot(log10(data.dose),SimFluoValues,'-', 'LineWidth', 2);
            errorbar(log10(data.dose),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI', 'FontSize',20)
            legend({'simulation'}, 'Location', 'southeast')
        end
    
    
            % From here on the steps are repeated as in the module (with
            % first repression coefficient) above
           
            %set the figure to be sceensize
            set(gcf, 'Position', get(0, 'Screensize'));
            %save figure
            saveas(gcf, ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_rep2_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.jpg']);
            %savefig(fig, '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_', num2str(draw))
            %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/SimulationFluoValues_draw_rep2_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.mat'], 'SimFluoValues');

end

for draw = 1:num_draws
    
     paraValues = table2array(rand_parameter(draw,:));
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

            NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues = simulate_DR_IPTG(para,data,ParaNames,model);
            DataMeans = data.data.means;
            DataStd = data.data.std;
            data.dose(1) = 1;
            
        %plot
        if FlPlot
            figure(3)
            subplot(4,3,draw);
            hold on;
            plot(log10(data.dose),SimFluoValues,'-', 'LineWidth', 2);
            errorbar(log10(data.dose),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI', 'FontSize',20)
            legend({'simulation'}, 'Location', 'southeast')
        end
    
    
            % From here on the steps are repeated as in the module (with
            % first repression coefficient) above
           
            %set the figure to be sceensize
            set(gcf, 'Position', get(0, 'Screensize'));
            %save figure
            saveas(gcf, ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_rep3_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.jpg']);
            %savefig(fig, '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_', num2str(draw))
            %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/SimulationFluoValues_draw_rep3_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.mat'], 'SimFluoValues');

end
            
        
        
