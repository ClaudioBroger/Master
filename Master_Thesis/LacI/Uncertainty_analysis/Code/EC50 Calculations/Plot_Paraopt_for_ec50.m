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

paraopt = array2table(paraopt);
paraopt.Properties.VariableNames = new_names;

paraopt.mu = 0.0077;
paraopt.kmaturation = 0.0173;
paraopt.dCit = 0;

model = sbioloadproject('LacImodel_uncertainty');
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
paraopt = paraopt(:,[string(ParaNames)]);

for k = 1:3
%% First repression coefficient

    paraValues = table2array(paraopt);
    paraValues = paraValues';
    
    objf  = [];
    objfn = {};
    ndata = 0;
    
    
    para = paraValues;
    dataPos = strcat(dataPath, "data.mat");
    data = load(dataPos);
    
    NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
    IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
    para(IdxToZero) = 0;
    
    SimFluoValues1_paraopt = simulate_DR_IPTG(para,data,ParaNames,model);
    DataMeans = data.data.means;
    DataStd = data.data.std;
    data.dose(1) = 1;
    
    %plot
    if FlPlot
        figure(1)
        hold on;
        plot(log10(data.dose),SimFluoValues1_paraopt,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-'),'Color',C(draw,:));
        errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
        xlabel('log IPTG (nM)', 'FontSize', 18)
        ylabel('mean Fluorescence','FontSize', 18)
        title('P4Lacn.2-cit + PAct1-LacI,optimal parameter 1st repression coefficient', 'FontSize',20)
    
    end
    hold off
    legend("show", 'Location', 'northeastoutside')

%% Second repression coefficient

    paraValues = table2array(paraopt);
    paraValues = paraValues';
    
    objf  = [];
    objfn = {};
    ndata = 0;
    
    
    para = paraValues;
    dataPos = strcat(dataPath, "data.mat");
    data = load(dataPos);
    
    NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep2', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
    IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
    para(IdxToZero) = 0;
    
    SimFluoValues2_paraopt = simulate_DR_IPTG(para,data,ParaNames,model);
    DataMeans = data.data.means;
    DataStd = data.data.std;
    data.dose(1) = 1;
    
    %plot
    if FlPlot
        figure(2)
        hold on;
        plot(log10(data.dose),SimFluoValues2_paraopt,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep2-'),'Color',C(draw,:));
        errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
        xlabel('log IPTG (nM)', 'FontSize', 18)
        ylabel('mean Fluorescence','FontSize', 18)
        title('P4Lacn.2-cit + PAct1-LacI, optimal parameter 2nd repression coefficient', 'FontSize',20)
    
    end
    hold off
    legend("show", 'Location', 'northeastoutside')

%% Third repression coefficient


    paraValues = table2array(paraopt);
    paraValues = paraValues';
    
    objf  = [];
    objfn = {};
    ndata = 0;
    
    
    para = paraValues;
    dataPos = strcat(dataPath, "data.mat");
    data = load(dataPos);
    
    NamestoZero = setdiff(ParaNames,{'kLacI','kCit', 'CitL', 'dLacI', 'LacIrep3', 'nLacI', 'KdLacI', 'mu', 'nMperUnit', 'kmaturation', 'indTime' });                            
    IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
    para(IdxToZero) = 0;
    
    SimFluoValues3_paraopt = simulate_DR_IPTG(para,data,ParaNames,model);
    DataMeans = data.data.means;
    DataStd = data.data.std;
    data.dose(1) = 1;
    
    %plot
    if FlPlot
        figure(3)
        hold on;
        plot(log10(data.dose),SimFluoValues3_paraopt,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep3-'),'Color',C(draw,:));
        errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
        xlabel('log IPTG (nM)', 'FontSize', 18)
        ylabel('mean Fluorescence','FontSize', 18)
        title('P4Lacn.2-cit + PAct1-LacI,optimal parameter 3rd repression coefficient', 'FontSize',20)
    
    end
    hold off
    legend("show", 'Location', 'northeastoutside')
end


SimFluoValues1_paraopt = array2table(SimFluoValues1_paraopt);
SimFluoValues2_paraopt = array2table(SimFluoValues2_paraopt);
SimFluoValues3_paraopt = array2table(SimFluoValues3_paraopt);

SimFluoValues1_paraopt(:,2) = array2table(data.dose);
SimFluoValues2_paraopt(:,2) = array2table(data.dose);
SimFluoValues3_paraopt(:,2) = array2table(data.dose);

SimFluoValues1_paraopt.Properties.VariableNames = ["SimFluoValues", "dose"];
SimFluoValues2_paraopt.Properties.VariableNames = ["SimFluoValues", "dose"];
SimFluoValues3_paraopt.Properties.VariableNames = ["SimFluoValues", "dose"];


