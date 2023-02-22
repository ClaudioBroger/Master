%set datapath to location of experimental data
dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
noise = false;
FlPlot = true;
for draw = 1:num_draws
    %find parameter names used in the model
    [ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
    %find parameter values from the model
    [paraValues] = cell2mat(get(model.mw_sbmod1.parameters, {'Value'}));
    %set the estimated parameter values to the random parameters
     paraValues(Settings.model.PIdx) = rand_parameter(:,draw);
    %TFlMode = FlMode;
    
     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;
    %Simulate every module three times; each time with another repression
    %coefficient -> 15 simulations per random parameter draw
        %% module 1: PAct1_LacI(WT)_tCyt1
                %set parameters (fixed and random drawn) to para
                para = paraValues;
                %load experimental data
                dataPos = strcat(dataPath, "data.mat");
                data = load(dataPos);
                
                %set parameters which are not used in this module to 0
                NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'LacI_rep_WT','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
                
                %simulate module using the fixed and randomly drawn
                %parameter values
                SimFluoValues(:,1) = simulate_DR_IPTG(para,data,ParaNames,model);
                DataMeans = data.data.means;
                DataStd = data.data.std;
                
                if noise           
                res = (SimFluoValues(:,1) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,1) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data.modulenames{1};
                ndata = ndata + numel(res);
                data.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,1);
                hold on;
                plot(log10(data.dose),SimFluoValues(:,1),'-', 'LineWidth', 2);
                %errorbar(log10(data.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence','FontSize', 18)
                title('P4Lacn.2-cit + PAct1-LacI', 'FontSize',20)
                legend({'simulation-WT'}, 'Location', 'southeast')
            end
    
    
            % From here on the steps are repeated as in the module (with
            % first repression coefficient) above
            %% module 1: PAct1_LacI(W220F)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data.mat");
                data_W220F = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'LacI_rep_W220F','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,2) = simulate_DR_IPTG(para,data_W220F,ParaNames,model);
                DataMeans = data_W220F.data.means;
                DataStd = data_W220F.data.std;
                
                if noise           
                res = (SimFluoValues(:,2) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,2) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,2);
                hold on;
                plot(log10(data_W220F.dose),SimFluoValues(:,2),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('PAct1-LacI(W220F)-tCyc1', 'FontSize', 20)
                legend({'simulation-W220F'}, 'Location', 'southeast')
            end
    
            %% module 1: PAct1_LacI(3mut)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data.mat");
                data_W220F = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'Silence_LacI_rep','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,3) = simulate_DR_IPTG(para,data_W220F,ParaNames,model);
                DataMeans = data_W220F.data.means;
                DataStd = data_W220F.data.std;
                
                if noise           
                res = (SimFluoValues(:,3) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,3) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,3);
                hold on;
                plot(log10(data_W220F.dose),SimFluoValues(:,3),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('PAct1-LacI(3mut)-tCyc1', 'FontSize', 20)
                legend({'simulation-3mut'}, 'Location', 'southeast')
            end
        
   
            %% module 2: PAct1_LacI(WT)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data_W220F.mat");
                data_W220F = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'LacI_rep_WT','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,4) = simulate_DR_IPTG(para,data_W220F,ParaNames,model);
                DataMeans = data_W220F.data.means;
                DataStd = data_W220F.data.std;
                
                if noise           
                res = (SimFluoValues(:,4) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,4) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,4);
                hold on;
                plot(log10(data_W220F.dose),SimFluoValues(:,4),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('PAct1-LacI(WT)-tCyc1', 'FontSize', 20)
                legend({'simulation-WT'}, 'Location', 'southeast')
            end
        
        

            %% module 2: PAct1_LacI(W220F)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data_W220F.mat");
                data_W220F = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'LacI_rep_W220F','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,5) = simulate_DR_IPTG(para,data_W220F,ParaNames,model);
                DataMeans = data_W220F.data.means;
                DataStd = data_W220F.data.std;
                
                if noise           
                res = (SimFluoValues(:,5) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,5) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,5);
                hold on;
                plot(log10(data_W220F.dose),SimFluoValues(:,5),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('PAct1-LacI(W220F)-tCyc1', 'FontSize', 20)
                legend({'simulation-W220F'}, 'Location', 'southeast')
            end
   
            %% module 2: PAct1_LacI(3mut)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data_W220F.mat");
                data_W220F = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'Silence_LacI_rep','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,6) = simulate_DR_IPTG(para,data_W220F,ParaNames,model);
                DataMeans = data_W220F.data.means;
                DataStd = data_W220F.data.std;
                
                if noise           
                res = (SimFluoValues(:,6) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,6) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,6);
                hold on;
                plot(log10(data_W220F.dose),SimFluoValues(:,6),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('PAct1-LacI(3mut)-tCyc1', 'FontSize', 20)
                legend({'simulation-3mut'}, 'Location', 'southeast')
            end
 
            %% module 3: P4Lacn.2_LacI(WT)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data_W220F_Q60G_T167A.mat");
                data_W220F_Q60G_T167A = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P_4Lacn_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'LacI_rep_WT','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,7) = simulate_DR_IPTG(para,data_W220F_Q60G_T167A,ParaNames,model);
                DataMeans = data_W220F_Q60G_T167A.data.means;
                DataStd = data_W220F_Q60G_T167A.data.std;
                
                if noise           
                res = (SimFluoValues(:,7) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,7) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F_Q60G_T167A.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F_Q60G_T167A.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,7);
                hold on;
                plot(log10(data_W220F_Q60G_T167A.dose),SimFluoValues(:,7),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F_Q60G_T167A.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P4Lacn.2-LacI(WT)-tCyc1', 'FontSize', 20)
                legend({'simulation-WT'}, 'Location', 'southeast')
            end
  
            %% module 3: P4Lacn.2_LacI(W220F)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data_W220F_Q60G_T167A.mat");
                data_W220F_Q60G_T167A = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P_4Lacn_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'LacI_rep_W220F','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,8) = simulate_DR_IPTG(para,data_W220F_Q60G_T167A,ParaNames,model);
                DataMeans = data_W220F_Q60G_T167A.data.means;
                DataStd = data_W220F_Q60G_T167A.data.std;
                
                if noise           
                res = (SimFluoValues(:,8) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,8) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F_Q60G_T167A.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F_Q60G_T167A.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,8);
                hold on;
                plot(log10(data_W220F_Q60G_T167A.dose),SimFluoValues(:,8),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F_Q60G_T167A.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P4Lacn.2-LacI(W220F)-tCyc1', 'FontSize', 20)
                legend({'simulation-W220F'}, 'Location', 'southeast')
            end
  
            %% module 3: P4Lacn.2_LacI(3mut)_tCyc1
                para = paraValues;
                dataPos = strcat(dataPath, "data_W220F_Q60G_T167A.mat");
                data_W220F_Q60G_T167A = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P_4Lacn_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'Silence_LacI_rep','LacI_rep_3mut','LacI_rep_3mut_P3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,9) = simulate_DR_IPTG(para,data_W220F_Q60G_T167A,ParaNames,model);
                DataMeans = data_W220F_Q60G_T167A.data.means;
                DataStd = data_W220F_Q60G_T167A.data.std;
                
                if noise           
                res = (SimFluoValues(:,9) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,9) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_W220F_Q60G_T167A.modulenames{1};
                ndata = ndata + numel(res);
                data_W220F_Q60G_T167A.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,9);
                hold on;
                plot(log10(data_W220F_Q60G_T167A.dose),SimFluoValues(:,9),'-', 'LineWidth', 2);
                %errorbar(log10(data_W220F_Q60G_T167A.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P4Lacn.2-LacI(3mut)-tCyc1', 'FontSize', 20)
                legend({'simulation-3mut'}, 'Location', 'southeast')
            end

            %% module 4: P3Lacn.5_LacI(WT)
                para = paraValues;
                dataPos = strcat(dataPath, "data_pt7.mat");
                data_pt7 = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P3_Lacn_5_cit','pt7_LacI', 'P3_Lacn_5_cit_L','LacI_rep_WT', 'dLacI_pt7', 'nLacI_P3', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,10) = simulate_DR_IPTG_3(para,data_pt7,ParaNames,model);
                DataMeans = data_pt7.data.means;
                DataStd = data_pt7.data.std;
                
                if noise           
                res = (SimFluoValues(:,10) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,10) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_pt7.modulenames{1};
                ndata = ndata + numel(res);
                data_pt7.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,10);
                hold on;
                plot(log10(data_pt7.dose),SimFluoValues(:,10),'-', 'LineWidth', 2);
                %errorbar(log10(data_pt7.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P3Lacn.5-LacI(WT)', 'FontSize', 20)
                legend({'simulation-WT'}, 'Location', 'southeast')
            end

            %% module 4: P3Lacn.5_LacI(W220F)
                para = paraValues;
                dataPos = strcat(dataPath, "data_pt7.mat");
                data_pt7 = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P3_Lacn_5_cit','pt7_LacI', 'P3_Lacn_5_cit_L','LacI_rep_W220F', 'dLacI_pt7', 'nLacI_P3', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,11) = simulate_DR_IPTG_3(para,data_pt7,ParaNames,model);
                DataMeans = data_pt7.data.means;
                DataStd = data_pt7.data.std;
                
                if noise           
                res = (SimFluoValues(:,11) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,11) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_pt7.modulenames{1};
                ndata = ndata + numel(res);
                data_pt7.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,11);
                hold on;
                plot(log10(data_pt7.dose),SimFluoValues(:,11),'-', 'LineWidth', 2);
                %errorbar(log10(data_pt7.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P3Lacn.5-LacI(W220F)', 'FontSize', 20)
                legend({'simulation-W220F'}, 'Location', 'southeast')
            end

            %% module 4: P3Lacn.5_LacI(3mut)
                para = paraValues;
                dataPos = strcat(dataPath, "data_pt7.mat");
                data_pt7 = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P3_Lacn_5_cit','pt7_LacI', 'P3_Lacn_5_cit_L','LacI_rep_3mut_P3','LacI_rep_3mut','Silence_LacI_rep', 'dLacI_pt7', 'nLacI_P3', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,12) = simulate_DR_IPTG_3(para,data_pt7,ParaNames,model);
                DataMeans = data_pt7.data.means;
                DataStd = data_pt7.data.std;
                
                if noise           
                res = (SimFluoValues(:,12) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,12) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_pt7.modulenames{1};
                ndata = ndata + numel(res);
                data_pt7.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,12);
                hold on;
                plot(log10(data_pt7.dose),SimFluoValues(:,12),'-', 'LineWidth', 2);
                %errorbar(log10(data_pt7.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P3Lacn.5-LacI(3mut)', 'FontSize', 20)
                legend({'simulation-3mut'}, 'Location', 'southeast')
            end
 
            %% module 5: P4Lacn.2_citrine_LacI(WT)
                para = paraValues;
                dataPos = strcat(dataPath, "data_pt7_5circuit.mat");
                data_pt7_5circuit = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P4Lacn_cit','pt7_LacI', 'P_4Lacn_LacI_L','LacI_rep_WT', 'dLacI_pt7', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,13) = simulate_DR_IPTG_3(para,data_pt7_5circuit,ParaNames,model);
                DataMeans = data_pt7_5circuit.data.means;
                DataStd = data_pt7_5circuit.data.std;
                
                if noise           
                res = (SimFluoValues(:,13) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,13) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_pt7_5circuit.modulenames{1};
                ndata = ndata + numel(res);
                data_pt7_5circuit.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,13);
                hold on;
                plot(log10(data_pt7_5circuit.dose),SimFluoValues(:,13),'-', 'LineWidth', 2);
                %errorbar(log10(data_pt7_5circuit.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P4Lacn.2-citrine-LacI(WT)', 'FontSize', 20)
                legend({'simulation-WT'}, 'Location', 'southeast')
            end

            %% module 5: P4Lacn.2_citrine_LacI(W220F)
                para = paraValues;
                dataPos = strcat(dataPath, "data_pt7_5circuit.mat");
                data_pt7_5circuit = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P4Lacn_cit','pt7_LacI', 'P_4Lacn_LacI_L','LacI_rep_W220F', 'dLacI_pt7', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,14) = simulate_DR_IPTG_3(para,data_pt7_5circuit,ParaNames,model);
                DataMeans = data_pt7_5circuit.data.means;
                DataStd = data_pt7_5circuit.data.std;
                
                if noise           
                res = (SimFluoValues(:,14) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,14) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_pt7_5circuit.modulenames{1};
                ndata = ndata + numel(res);
                data_pt7_5circuit.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,14);
                hold on;
                plot(log10(data_pt7_5circuit.dose),SimFluoValues(:,14),'-', 'LineWidth', 2);
                %errorbar(log10(data_pt7_5circuit.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P4Lacn.2-citrine-LacI(W220F)', 'FontSize', 20)
                legend({'simulation-W220F'}, 'Location', 'southeast')
            end

            %% module 5: P4Lacn.2_citrine_LacI(3mut)
                para = paraValues;
                dataPos = strcat(dataPath, "data_pt7_5circuit.mat");
                data_pt7_5circuit = load(dataPos);
    
                NamestoZero = setdiff(ParaNames,{'P4Lacn_cit','pt7_LacI', 'P_4Lacn_LacI_L','LacI_rep_3mut','LacI_rep_3mut_P3','Silence_LacI_rep', 'dLacI_pt7', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
                IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
                para(IdxToZero) = 0;
    
                SimFluoValues(:,15) = simulate_DR_IPTG_3(para,data_pt7_5circuit,ParaNames,model);
                DataMeans = data_pt7_5circuit.data.means;
                DataStd = data_pt7_5circuit.data.std;
                
                if noise           
                res = (SimFluoValues(:,15) + normrnd(0,DataStd) - DataMeans)./DataStd;
                else
                res = (SimFluoValues(:,15) - DataMeans)./DataStd;
                end
     
                objf(end+1) = sum(res.^2);
                objfn{end+1} = data_pt7_5circuit.modulenames{1};
                ndata = ndata + numel(res);
                data_pt7_5circuit.dose(1) = 1;
            %plot
            if FlPlot
                figure(draw);
                subplot(5,3,15);
                hold on;
                plot(log10(data_pt7_5circuit.dose),SimFluoValues(:,15),'-', 'LineWidth', 2);
                %errorbar(log10(data_pt7_5circuit.dose),DataMeans,DataStd,'o');
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence', 'FontSize', 18)
                title('P4Lacn.2-citrine-LacI(3mut)', 'FontSize', 20)
                legend({'simulation-3mut'}, 'Location', 'southeast')
            end
            %set the figure to be sceensize
            set(gcf, 'Position', get(0, 'Screensize'));
            %save figure
            saveas(gcf, ['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.jpg']);
            %savefig(fig, '/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Figure_draw_', num2str(draw))
            %save simulated fluoresence values
            save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/SimulationFluoValues_draw_' num2str(draw),'_',datestr(now, 'dd-mmm-yyyy'), '.mat'], 'SimFluoValues');
            delete(gcf)
end
        
        
