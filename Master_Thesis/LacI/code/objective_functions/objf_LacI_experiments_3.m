function [objfc, ndata, objf, objfn] = objf_LacI_experiments_3(paraest, model, dataPath, FlPlot, FlMode, PIdx, noise)
% objf_TetR_experiments: simulates experiments and compares simulation to
% data.
% 
% Output Arguments
% =================
% - objfc: total cost
% - ndata: number of data points
% - objf: vector of costs for each separate experiment
% - objfn: labels of the objf vectors

    % INITIALIZE ARGUMENTS %
    %%%%%%%%%%%%%%%%%%%%%%%%

    [ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
    [paraValues] = cell2mat(get(model.mw_sbmod1.parameters, {'Value'}));
     paraValues(PIdx) = paraest;
     TFlMode = FlMode;
 
    % INITIALIZE OUTPUTS %
    %%%%%%%%%%%%%%%%%%%%%%
    objf  = [];
    objfn = {};
    ndata = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~isempty(TFlMode) %true when TflMode is not empty
    TFl = TFlMode(1);
    switch TFl              
        case{1}
        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            dataPos = strcat(dataPath, "data.mat");
            data = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'dCit', 'LacI_rep_WT','LacI_rep_3mut', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues = simulate_DR_IPTG(para,data.data,ParaNames,model);
            DataMeans = data.data.means;
            DataStd = data.data.std;
            
            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end
 
            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);
        %plot
        if FlPlot
            figure(1)
            subplot(2,2,1);
            hold on;
            plot(log10(data.dose),SimFluoValues,'-');
            errorbar(log10(data.dose),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)')
            ylabel('mean Fluorescence')
            title('P4Lacn.2_cit + PAct1_LacI')
            legend({'simulation' 'data'})
        end
        case{2}
        %% module 2: PAct1_LacI(W220F)_tCyc1
            para = paraValues;
            dataPos = strcat(dataPath, "data_W220F.mat");
            data_W220F = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'PAct1_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'dCit', 'LacI_rep_W220F','LacI_rep_3mut', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues = simulate_DR_IPTG(para,data_W220F.data_W220F,ParaNames,model);
            DataMeans = data_W220F.data_W220F.means;
            DataStd = data_W220F.data_W220F.std;
            
            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end
 
            objf(end+1) = sum(res.^2);
            objfn{end+1} = data_W220F.modulenames_W220F{1};
            ndata = ndata + numel(res);
        %plot
        if FlPlot
            figure(1)
            subplot(2,2,2);
            hold on;
            plot(log10(data_W220F.dose_W220F),SimFluoValues,'-');
            errorbar(log10(data_W220F.dose_W220F),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)')
            ylabel('mean Fluorescence')
            title('PAct1_LacI(W220F)_tCyc1')
            legend({'simulation_W220F' 'data_W220F'})
        end
        case{3}
        %% module 3: P4Lacn.2_LacI(W220F,Q60G,T167A)_tCyc1
            para = paraValues;
            dataPos = strcat(dataPath, "data_W220F_Q60G_T167A.mat");
            data_W220F_Q60G_T167A = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'P_4Lacn_LacI','P4Lacn_cit', 'P_4Lacn_LacI_L', 'dLacI', 'dCit', 'Silence_LacI_rep','LacI_rep_3mut', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues = simulate_DR_IPTG(para,data_W220F_Q60G_T167A.data_W220F_Q60G_T167A,ParaNames,model);
            DataMeans = data_W220F_Q60G_T167A.data_W220F_Q60G_T167A.means;
            DataStd = data_W220F_Q60G_T167A.data_W220F_Q60G_T167A.std;
            
            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end
 
            objf(end+1) = sum(res.^2);
            objfn{end+1} = data_W220F_Q60G_T167A.modulenames_W220F_Q60G_T167A{1};
            ndata = ndata + numel(res);
        %plot
        if FlPlot
            figure(1)
            subplot(2,2,3);
            hold on;
            plot(log10(data_W220F_Q60G_T167A.dose_W220F_Q60G_T167A),SimFluoValues,'-');
            errorbar(log10(data_W220F_Q60G_T167A.dose_W220F_Q60G_T167A),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)')
            ylabel('mean Fluorescence')
            title('P4Lacn.2_LacI(W220F,Q60G, T167A)_tCyc1')
            legend({'simulation_W220F_Q60G_T167A' 'data_W220F_Q60G_T167A'})
        end
        case{4}
        %% module 4: P3Lacn.5_LacI(W220F,Q60G,T167A)
            para = paraValues;
            dataPos = strcat(dataPath, "data_pt7.mat");
            data_pt7 = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'P3_Lacn_5_cit','pt7_LacI', 'P3_Lacn_5_cit_L','LacI_rep_3mut','Silence_LacI_rep', 'dLacI_pt7', 'dCit', 'nLacI_P3', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues = simulate_DR_IPTG_3(para,data_pt7.data_pt7,ParaNames,model);
            DataMeans = data_pt7.data_pt7.means;
            DataStd = data_pt7.data_pt7.std;
            
            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end
 
            objf(end+1) = sum(res.^2);
            objfn{end+1} = data_pt7.modulenames_pt7{1};
            ndata = ndata + numel(res);
        %plot
        if FlPlot
            figure(1)
            subplot(2,2,4);
            hold on;
            plot(log10(data_pt7.dose_pt7),SimFluoValues,'-');
            errorbar(log10(data_pt7.dose_pt7),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)')
            ylabel('mean Fluorescence')
            title('P3Lacn.5_LacI(W220F,Q60G,T167A)')
            legend({'simulation_pt7' 'data_pt7'})
        end
    end 
    
    TFlMode = setdiff(TFlMode,TFl); %removes the mode that has just been dealt with
end

objfc = sum(objf); %total cost

end