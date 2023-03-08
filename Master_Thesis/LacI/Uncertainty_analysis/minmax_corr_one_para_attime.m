dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');

random = 1;
around_paraopt = 0;

dataPos = strcat(dataPath, "data.mat");
data = load(dataPos);

for k= 1:num_draws
        figure(1)
        subplot(4,5,k)
        plot(data.data.means(:), SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k))
        c = corrcoef(data.data.means(:), SimFluoValues1combined.SimFluoValues(SimFluoValues1combined.Draw == k));
        corr(k,1) = c(1,2);
end

index_min = find(corr == min(corr));
index_max = find(corr == max(corr));

if random
    for k = 1:width(rand_parameter)
        diff(k) = table2array(rand_parameter(index_min, k)) - table2array(rand_parameter(index_max, k));
    end
    rand_parameter = table2array(rand_parameter);
    for k = 1:width(rand_parameter)
        max_range(k) = max(rand_parameter(:,k));
        min_range(k) = min(rand_parameter(:,k));
    end
    
    
    
    
    objf  = [];
    objfn = {};
    ndata = 0;

    
    C = linspecer(width(rand_parameter));
    for variable = 1:width(rand_parameter)
        %%Min parameter 
        parametervalues = rand_parameter(index_max,:);
        parametervalues(variable) = min_range(variable);
    
    
        para = parametervalues';
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
            figure(2)
            hold on;
            plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-min-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep1 1 parameter set to minimum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    
        %%max parameter
        parametervalues = rand_parameter(index_max,:);
        parametervalues(variable) = max_range(variable);
    
    
        para = parametervalues';
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
            figure(3)
            hold on;
            plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-max-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep1 1 parameter set to maximum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    end
end

if around_paraopt
    for k = 1:width(paraopt_range)
        diff(k) = table2array(paraopt_range(index_min, k)) - table2array(paraopt_range(index_max, k));
    end
    paraopt_range = table2array(paraopt_range);
    for k = 1:width(paraopt_range)
        max_range(k) = max(paraopt_range(:,k));
        min_range(k) = min(paraopt_range(:,k));
    end
    
    
    
    
    objf  = [];
    objfn = {};
    ndata = 0;

    
    C = linspecer(width(paraopt_range));
    for variable = 1:width(paraopt_range)
        %%Min parameter 
        parametervalues = paraopt_range(index_max,:);
        parametervalues(variable) = min_range(variable);
    
    
        para = parametervalues';
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
            figure(4)
            hold on;
            plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-min-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 1 parameter set to minimum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    
        %%max parameter
        parametervalues = paraopt_range(index_max,:);
        parametervalues(variable) = max_range(variable);
    
    
        para = parametervalues';
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
            figure(5)
            hold on;
            plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-max-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 1 parameter set to maximum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    end
end




for k= 1:num_draws

        subplot(4,5,k)
        figure(6)
        plot(data.data.means(:), SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k))
        c = corrcoef(data.data.means(:), SimFluoValues2combined.SimFluoValues(SimFluoValues2combined.Draw == k));
        corr2(k,1) = c(1,2);
end

index_min = find(corr2 == min(corr2));
index_max = find(corr2 == max(corr2));

if random
    for k = 1:width(rand_parameter)
        diff(k) = rand_parameter(index_min, k) - rand_parameter(index_max, k);
    end

    for k = 1:width(rand_parameter)
        max_range(k) = max(rand_parameter(:,k));
        min_range(k) = min(rand_parameter(:,k));
    end
    
    
    
    
    objf  = [];
    objfn = {};
    ndata = 0;

    
    C = linspecer(width(rand_parameter));
    for variable = 1:width(rand_parameter)
        %%Min parameter 
        parametervalues = rand_parameter(index_max,:);
        parametervalues(variable) = min_range(variable);
    
    
        para = parametervalues';
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
            figure(7)
            hold on;
            plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-min-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep2 1 parameter set to minimum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    
        %%max parameter
        parametervalues = rand_parameter(index_max,:);
        parametervalues(variable) = max_range(variable);
    
    
        para = parametervalues';
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
            figure(8)
            hold on;
            plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-max-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep2 1 parameter set to maximum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    end
end

if around_paraopt
    for k = 1:width(paraopt_range)
        diff(k) = paraopt_range(index_min, k) - paraopt_range(index_max, k);
    end

    for k = 1:width(paraopt_range)
        max_range(k) = max(paraopt_range(:,k));
        min_range(k) = min(paraopt_range(:,k));
    end
    
    
    
    
    objf  = [];
    objfn = {};
    ndata = 0;

    
    C = linspecer(width(paraopt_range));
    for variable = 1:width(paraopt_range)
        %%Min parameter 
        parametervalues = paraopt_range(index_max,:);
        parametervalues(variable) = min_range(variable);
    
    
        para = parametervalues';
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
            figure(9)
            hold on;
            plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-min-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep2 1 parameter set to minimum at the time', 'FontSize',20)
            
        end
        hold off
        legend("show", 'Location', 'northeastoutside')
    
        %%max parameter
        parametervalues = paraopt_range(index_max,:);
        parametervalues(variable) = max_range(variable);
    
    
        para = parametervalues';
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
            figure(10)
            hold on;
            plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-max-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep2 1 parameter set to maximum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    end
end


%% 3rd LacI repression coefficient
for k= 1:num_draws

        subplot(4,5,k)
        figure(11)
        plot(data.data.means(:), SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k))
        c = corrcoef(data.data.means(:), SimFluoValues3combined.SimFluoValues(SimFluoValues3combined.Draw == k));
        corr3(k,1) = c(1,2);
end

index_min = find(corr3 == min(corr3));
index_max = find(corr3 == max(corr3));

if random
    for k = 1:width(rand_parameter)
        diff(k) = rand_parameter(index_min, k) - rand_parameter(index_max, k);
    end

    for k = 1:width(rand_parameter)
        max_range(k) = max(rand_parameter(:,k));
        min_range(k) = min(rand_parameter(:,k));
    end
    
    
    
    
    objf  = [];
    objfn = {};
    ndata = 0;

    
    C = linspecer(width(rand_parameter));
    for variable = 1:width(rand_parameter)
        %%Min parameter 
        parametervalues = rand_parameter(index_max,:);
        parametervalues(variable) = min_range(variable);
    
    
        para = parametervalues';
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
            figure(12)
            hold on;
            plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-min-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep3 1 parameter set to minimum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    
        %%max parameter
        parametervalues = rand_parameter(index_max,:);
        parametervalues(variable) = max_range(variable);
    
    
        para = parametervalues';
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
            figure(13)
            hold on;
            plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-max-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep3 1 parameter set to maximum at the time', 'FontSize',20)
            
        end
        legend("show", 'Location', 'northeastoutside')
    end
end

if around_paraopt
    for k = 1:width(paraopt_range)
        diff(k) = paraopt_range(index_min, k) - paraopt_range(index_max, k);
    end

    for k = 1:width(paraopt_range)
        max_range(k) = max(paraopt_range(:,k));
        min_range(k) = min(paraopt_range(:,k));
    end
    
    
    
    
    objf  = [];
    objfn = {};
    ndata = 0;

    
    C = linspecer(width(paraopt_range));
    for variable = 1:width(paraopt_range)
        %%Min parameter 
        parametervalues = paraopt_range(index_max,:);
        parametervalues(variable) = min_range(variable);
    
    
        para = parametervalues';
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
            figure(14)
            hold on;
            plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-min-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep3 1 parameter set to minimum at the time', 'FontSize',20)
            
        end
        hold off
        legend("show", 'Location', 'northeastoutside')
    
        %%max parameter
        parametervalues = paraopt_range(index_max,:);
        parametervalues(variable) = max_range(variable);
    
    
        para = parametervalues';
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
            figure(15)
            hold on;
            plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-max-parameter-',string(ParaNames(variable))), 'Color', C(variable,:));
            errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI,LacIrep3 1 parameter set to maximum at the time', 'FontSize',20)
            
        end
        hold off
        legend("show", 'Location', 'northeastoutside')
    end
end
