function [objfc, ndata, objf, objfn] = objfEvolvexp(paraest, model, dataPath, FlPlot, FlMode, PIdx, noise)
    % objf: simulates experiments and compares simulation to
    % data. Added minstdv
    % % actual experiment 
    % 
    % Output Arguments
    % =================
    % - objfc: total cost
    % - ndata: number of data points
    % - objf: vector of costs for each separate experiment
    % - objfn: labels of the objf vectors
    


    % INITIALIZE ARGUMENTS %
    %%%%%%%%%%%%%%%%%%%%%%%%

    [ParaNames, paraValues]   = IQMparameters(model);
    paraValues(PIdx) = paraest;
    TFlMode = FlMode;
    
    
    % INITIALIZE OUTPUTS %
    %%%%%%%%%%%%%%%%%%%%%%
    objf  = [];
    objfn = {};
    ndata = 0;

%     if FlPlot
%         figure;
%         set(gcf,'color','w');
%     end
%     
    min_stdv = 1.1163;
    idxGRmax = find(ismember(ParaNames, 'growthMAX'));
    
    % EVALUATE %
    %%%%%%%%%%%%
    while ~isempty(TFlMode)
        TFl = TFlMode(1);
        switch TFl  
            case{1}
            %% module 1: pact1_TetR + p7tet1_Cit (Simple Repression)
            para = paraValues;
            dataPos = strcat(dataPath, "dataSRtet.mat");
            data = load(dataPos);
            data.data.time = 420;
            NamestoZero = setdiff(ParaNames,{'kactTetR','k7tetCit', 'kL7tet', 'dTetR', 'dCit', 'thetaTetR', 'nTetR', 'f', 'nTup1', 'KdTetR', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime' , 'thetaTup1', 'thetaLacI3mut_2', 'a', 'thetaLacI3mut', 'p1', 'p2'});                
            IdxToZero = find(ismember(ParaNames, NamestoZero));           
            para(IdxToZero) = 0;                  

%             
            %simulate the fluo values using the timepoints from the data
            SimFluoValues = simulateDRaTc(para,data.data,ParaNames,model);

            %calculate the cost of the simulation
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end

            % evaluate outputs
            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            %plot
            if FlPlot
                subplot(3,4,1);
                hold on;
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);               
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Fluorescence (a.u.)')
                %title(data.modulenames{1})
                title('Simple Repression TetR')
                legend({ 'data' 'simulation'},'Location','northwest')
                hold off
            end
            case{2}
            %% module 2: p7tet1_TetR + p7tet1_Cit (Auto Repression)
            para = paraValues;
            dataPos = strcat(dataPath, "dataARtet.mat");
            data = load(dataPos);
            data.data.time = 420;
            NamestoZero = setdiff(ParaNames,{'k7tetTetR','k7tetCit', 'kL7tet', 'dTetR', 'dCit', 'thetaTetR', 'nTetR','f', 'nTup1', 'KdTetR', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime' , 'thetaTup1', 'thetaLacI3mut_2', 'a', 'thetaLacI3mut', 'p1', 'p2'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues = simulateDRaTc(para,data.data,ParaNames,model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            
            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,2);
                hold on;
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
               
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Fluorescence (a.u.)');
                %title(data.modulenames{1})
                title('Autorepression TetR')
                %legend({'simulation' 'data'},'Location','southeast');
            end
            case{3}
            %% module 3: pRNR2_TetRTup1 + p7tet1_Cit
            para = paraValues;
            dataPos = strcat(dataPath, "dataSRtup1.mat");
            data = load(dataPos);
            data.data.time = 360;
            NamestoZero = setdiff(ParaNames,{'krnrTup1','k7tetCit', 'kL7tet', 'dTup1', 'dCit', 'thetaTup1',  'nTup1', 'g', 'nTetR' 'KdTetR', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime' ,'thetaTetR','thetaLacI3mut_2', 'a', 'thetaLacI3mut', 'p1', 'p2'});                
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues = simulateDRaTc(para,data.data,ParaNames,model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            
            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,3);
                hold on;
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Fluorescence (a.u.)')
                %title(data.modulenames{1})
                title('Simple Repression TetR-Tup1')
                %legend({'simulation' 'data'},'Location','southeast')
                hold off
            end
            case{4}
            %% module 4: WTC p7tet1_TetR + pRNR2_TetRTup1 + p7tet1_Cit
            para = paraValues;
            dataPos = strcat(dataPath, "dataWTCcit.mat");
            data = load(dataPos);

            NamestoZero = setdiff(ParaNames,{'krnrTup1', 'k7tetTetR','k7tetCit', 'kL7tet', 'dTetR', 'dCit', 'dTup1', 'f', 'g', 'thetaTetR', 'thetaTup1', 'nTetR', 'nTup1', 'KdTetR', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime', 'thetaLacI3mut_2', 'a', 'thetaLacI3mut', 'p1', 'p2'});                
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            data.data.time = 360;
            SimFluoValues = simulateDRaTc(para,data.data,ParaNames,model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);

            if noise           
            res = (SimFluoValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues - DataMeans)./DataStd;
            end

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,4);
                hold on;
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
               
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Fluorescence (a.u.)')
                %title(data.modulenames{1})
                title('WTC')
                %legend({'simulation' 'data'},'Location','northwest')
                hold off
            end
            case{5}
            %% module 5: WTC p7tet1_TetR + pRNR2_TetRTup1 + p7tet1_PMA1
            para = paraValues;
            dataPos = strcat(dataPath, "dataWTCpma1.mat");
            data = load(dataPos);
            %data.data.time = 1440;
            NamestoZero = setdiff(ParaNames,{'krnrTup1', 'k7tetTetR','k7tetPMA1', 'kL7tet', 'dTetR', 'dTup1', 'dPMA1', 'thetaTetR', 'thetaTup1', 'nTetR', 'nTup1', 'f', 'g', 'KdTetR', 'p1', 'p2', 'growthMIN', 'growthMAX' 'nMperUnit', 'indTime', 'thetaLacI3mut_2', 'a', 'thetaLacI3mut'});        
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimGrowthValues  = simulateDRgrowth(para,data.data,ParaNames, model);
            DataMeans = data.data.means.';
            DataStd = data.data.std.';
            if noise           
            res = (SimGrowthValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimGrowthValues - DataMeans)./DataStd;
            end
            

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,5);
                hold on;
                %plot(log10(data.data.dose),DataMeans,'-','Linewidth',5);
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimGrowthValues,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Growth rate (min^{-1})')
                %title(data.modulenames{1})
                title('WTC PMA1')
                %legend({'simulation' 'data'},'Location','northwest')
                hold off
            end
             case{6}
            %% module 6: P4Lacn.2_cit + PAct1_LacI
            para = paraValues;
            dataPos = strcat(dataPath, "dataSRlacI.mat");
            data = load(dataPos);
            data.data.time = 1200;
            NamestoZero = setdiff(ParaNames,{'kactLacI', 'kl4lac', 'dLacI', 'nLacI', 'KdLacI', 'thetaLacI', 'k4lacCit', 'dCit', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime', 'thetaLacI3mut_2', 'b', 'thetaLacI3mut', 'thetaTup1', 'thetaTetR', 'p1', 'p2', 'nTetR', 'nTup1'});        
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues = simulateDRiptg(para,data.data,ParaNames, model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            if noise           
            res = (SimFluoValues  + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues  - DataMeans)./DataStd;
            end
            

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,6);
                hold on;
                %plot(log10(data.data.dose),DataMeans,'-','Linewidth',5);
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues ,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Growth rate (min^{-1})')
                %title(data.modulenames{1})
                title('SR LacI')
                %legend({'simulation' 'data'},'Location','northwest')
                hold off
            end
        case{7}
            %% module 7: P4Lacn.2_cit + PAct1_LacI_W22OF
            para = paraValues;
            dataPos = strcat(dataPath, "dataSRlacI_mut.mat");
            data = load(dataPos);
            data.data.time = 1200;
            NamestoZero = setdiff(ParaNames,{'kactLacI', 'kl4lac', 'dLacI', 'nLacI', 'KdLacI', 'thetaLacImut', 'k4lacCit', 'dCit', 'growthFix', 'nMperUnit', 'kmaturation', 'b', 'indTime', 'thetaLacI3mut_2', 'thetaLacI3mut', 'thetaTup1', 'thetaTetR', 'p1', 'p2', 'nTetR', 'nTup1', });        
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues = simulateDRiptg(para,data.data,ParaNames, model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            if noise           
            res = (SimFluoValues  + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues  - DataMeans)./DataStd;
            end
            

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,7);
                hold on;
                %plot(log10(data.data.dose),DataMeans,'-','Linewidth',5);
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues ,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Growth rate (min^{-1})')
                %title(data.modulenames{1})
                title('SR LacI 1x mut')
                %legend({'simulation' 'data'},'Location','northwest')
                hold off
            end
                 case{8}
            %% module 8: P4Lacn.2_cit + P4Lacn.2_LacI_3xmut
            para = paraValues;
            dataPos = strcat(dataPath, "dataARlacI_mut3.mat");
            data = load(dataPos);
            data.data.time = 1200;
            NamestoZero = setdiff(ParaNames,{'k4lacLacI', 'kl4lac', 'dLacI', 'nLacI', 'KdLacI', 'thetaLacI3mut', 'a', 'k4lacCit', 'dCit', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime', 'thetaLacI3mut_2', 'b', 'thetaTup1', 'thetaTetR', 'p1', 'p2', 'nTetR', 'nTup1', });        
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues   = simulateDRiptg(para,data.data,ParaNames, model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            if noise           
            res = (SimFluoValues  + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues  - DataMeans)./DataStd;
            end
            

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,8);
                hold on;
                %plot(log10(data.data.dose),DataMeans,'-','Linewidth',5);
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues ,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Growth rate (min^{-1})')
                %title(data.modulenames{1})
                title('AR LacI 3 mut')
                %legend({'simulation' 'data'},'Location','northwest')
                hold off
            end
         case{9}
            %% module 9: Pact_Cit_LacI + P3Lacn.2_LacI_3xmut
            para = paraValues;
            dataPos = strcat(dataPath, "dataSRlacI_mut3.mat");
            data = load(dataPos);
            data.data.time = 1200;
            NamestoZero = setdiff(ParaNames,{'kactLacICit', 'kl3lac', 'dLacICit', 'nLacI2', 'KdLacI', 'k3lacCit', 'dCit', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime', 'thetaLacI3mut_2', 'thetaLacI3mut', 'a', 'b', 'thetaTup1', 'thetaTetR', 'p1', 'p2', 'nTetR', 'nTup1' });        
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues  = simulateDRiptg_dblCit(para,data.data,ParaNames, model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            if noise           
            res = (SimFluoValues  + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues  - DataMeans)./DataStd;
            end
            

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,9);
                hold on;
                %plot(log10(data.data.dose),DataMeans,'-','Linewidth',5);
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues ,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Growth rate (min^{-1})')
                %title(data.modulenames{1})
                title('AR LacI_Cit 3 mut, k3lac')
                %legend({'simulation' 'data'},'Location','northwest')
                hold off
            end  
        case{10}
            %% module 10: Pact_Cit_LacI + P4Lacn.2_LacI_3xmut
            para = paraValues;
            dataPos = strcat(dataPath, "dataSRlacI_mut4.mat");
            data = load(dataPos);
            data.data.time = 1200;
            NamestoZero = setdiff(ParaNames,{'kactLacICit', 'kl4lac', 'dLacICit', 'nLacI', 'KdLacI', 'k4lacCit', 'thetaLacI3mut', 'a', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime', 'thetaLacI3mut_2', 'b', 'thetaTup1', 'thetaTetR', 'p1', 'p2', 'nTetR', 'nTup1', });        
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

           SimFluoValues  = simulateDRiptg_dblCit(para,data.data,ParaNames, model);
            DataMeans = data.data.means;
            %DataStd = data.data.std;
            DataStd = max(min_stdv * ones(size(DataMeans)), sqrt(data.data.std + data.data.empty_stdv).^2);
            if noise           
            res = (SimFluoValues  + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimFluoValues  - DataMeans)./DataStd;
            end
            

            objf(end+1) = sum(res.^2);
            objfn{end+1} = data.modulenames{1};
            ndata = ndata + numel(res);

            if FlPlot
                subplot(3,4,10);
                hold on;
                %plot(log10(data.data.dose),DataMeans,'-','Linewidth',5);
                errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimFluoValues ,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)')
                ylabel('mean Growth rate (min^{-1})')
                %title(data.modulenames{1})
                title('AR LacI-Cit 3 mut, k4lac')
                %legend({'simulation' 'data'},'Location','northwest')
                hold off
            end   
            
     case{11}        
              %% module 11: pact1_cit-Lac + p3lacn5-tup1-adh1-noNLS (2x) + p7tet1.PMA1 
              
            dataPos = strcat(dataPath, "dataDDRLacIPMA1.mat");
            data = load(dataPos);
            para = paraValues;
            data.data.dose = (data.data.dose/468.22)*1000; 
            
            NamestoZero = setdiff(ParaNames,{'k3lacTup1', 'k7tetPMA1', 'kactLacICit' , 'kl3lac', 'dLacICit', 'dPMA1',  'dTup1',  'dTup1adh1', 'thetaLacI3mut', 'thetaLacI3mut_2', 'a', 'g', 'nTup1', 'nLacI2', 'KdTetR', 'KdLacI', 'nMperUnit', 'indTime' , 'thetaTup1', 'thetaLacI_2', 'p1', 'p2', 'thetaTetR', 'nTetR', 'growthMIN', 'growthMAX'});                
            IdxToZero = find(ismember(ParaNames, NamestoZero));           
            para(IdxToZero) = 0; 
            para(idxGRmax) = 0.0085;  
                  
            IdxTup = find(ismember(ParaNames, 'k3lacTup1'));
            para(IdxTup) = 2*para(IdxTup);                  
            IdxTupLeak = find(ismember(ParaNames, 'kl3lac'));
            para(IdxTupLeak) = 2*para(IdxTupLeak);
            
            %simulate the fluo values using the timepoints from the data
            SimGrowthValues =  simulateDRgrowth(para, data.data, ParaNames,model);
            DataMeans = data.data.means.';
            DataStd = data.data.std.';
            %data.data.dose(1) = data.data.dose(1) + 0.00001;
            if noise           
            res = (SimGrowthValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimGrowthValues - DataMeans)./(DataStd);
            end

            objf(end+1) = sum(sum(res.^2));
            objfn{end+1} = 'Test Growth';
            ndata = ndata + numel(res); 

            if FlPlot
                subplot(3,4,11);               
                hold on;
               errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimGrowthValues,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);             
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)', 'fontweight','bold')
                ylabel('mean Growth Rate (a.u.)', 'fontweight','bold')
                title('DDR LacI-Tup1(2x)-PMA1')
                %legend({'simulation' 'data'},'Location','southeast')
                hold off
            end 
     case{12}        
              %% module 12: pact1_cit-Lac-mut3x + p3lacn5-tup1-adh1-noNLS (2x) + p7tet1.PMA1 
              
            dataPos = strcat(dataPath, "dataDDRLacIMutPMA1.mat");
            data = load(dataPos);
            para = paraValues;
            data.data.dose = (data.data.dose/468.22)*1000; 
            
            NamestoZero = setdiff(ParaNames,{'k3lacTup1', 'k7tetPMA1', 'kactLacICit' , 'kl3lac', 'thetaLacI3mut_2', 'dLacICit', 'dPMA1',  'dTup1',  'dTup1adh1', 'a', 'b', 'g', 'nTup1', 'nLacI2', 'KdTetR', 'KdLacI', 'nMperUnit', 'thetaLacI3mut', 'thetaTetR', 'nTetR', 'indTime' , 'thetaTup1', 'p1', 'p2', 'growthMIN', 'growthMAX'});                
            IdxToZero = find(ismember(ParaNames, NamestoZero));           
            para(IdxToZero) = 0; 
            para(idxGRmax) = 0.0081; 
                  
            IdxTup = find(ismember(ParaNames, 'k3lacTup1'));
            para(IdxTup) = 2*para(IdxTup);                  
            IdxTupLeak = find(ismember(ParaNames, 'kl3lac'));
            para(IdxTupLeak) = 2*para(IdxTupLeak);
            
            %simulate the fluo values using the timepoints from the data
            SimGrowthValues =  simulateDRgrowth(para, data.data, ParaNames,model);
            DataMeans = data.data.means.';
            DataStd = data.data.std.';
           
            if noise           
            res = (SimGrowthValues + normrnd(0,DataStd) - DataMeans)./DataStd;
            else
            res = (SimGrowthValues - DataMeans)./(DataStd);
            end

            objf(end+1) = sum(sum(res.^2));
            objfn{end+1} = 'Test Growth';
            ndata = ndata + numel(res); 

            %data.data.dose(1) = data.data.dose(1) + 0.00001;
            if FlPlot
                subplot(3,4,12);               
                hold on;
              errorbar(log10(data.data.dose),DataMeans,DataStd,'o','Color',[0, 0.4470, 0.7410],'Linewidth',1.5);
                plot(log10(data.data.dose),SimGrowthValues,'-','Linewidth',2, 'Color', [0.9100 0.4100 0.1700]);
                
                set(gca,'Fontsize',7);
                xlabel('log_{10} aTc (nM)', 'fontweight','bold')
                ylabel('mean Growth Rate (a.u.)', 'fontweight','bold')
                title('DDR LacI (mut)-Tup1(2x)-PMA1')
                %legend({'simulation' 'data'},'Location','southeast')
                hold off
            end 
        end    

        TFlMode = setdiff(TFlMode,TFl);
    end

    objfc = sum(objf);
    %objfc = -(1/2)*sum(objf);
end