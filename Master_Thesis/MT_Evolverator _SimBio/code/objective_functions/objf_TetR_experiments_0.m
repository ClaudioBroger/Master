function [objfc, ndata, objf, objfn] = objf_TetR_experiments_0(paraest, model, dataPath, FlPlot, FlMode, PIdx, noise)
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
            dataPos = strcat(dataPath, "dataARtet.mat");
            data = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'k7tetTetR','k7tetCit', 'kL7tet', 'dTetR', 'dCit', 'thetaTetR', 'nTetR', 'KdTetR', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues = simulate_DR_aTc(para,data.data,ParaNames,model);
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
            hold on;
            plot(log10(data.data.dose),SimFluoValues,'-');
            errorbar(log10(data.data.dose),DataMeans,DataStd,'o');
            xlabel('log aTc (nM)')
            ylabel('mean Fluorescence')
            title('ptetTetR-ptetCitrine')
            legend({'simulation' 'data'})
        end
    end 
    
    TFlMode = setdiff(TFlMode,TFl); %removes the mode that has just been dealt with
end

objfc = sum(objf); %total cost

end