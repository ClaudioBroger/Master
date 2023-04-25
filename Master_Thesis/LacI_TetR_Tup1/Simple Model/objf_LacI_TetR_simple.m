function [objfc, ndata, objf, objfn, SimFluoValues1, SimFluoValues2, SimFluoValues3] = objf_LacI_TetR(paraest, model, PIdx, noise, paraValues, ParaNames, data)
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

   
 
    % INITIALIZE OUTPUTS %
    %%%%%%%%%%%%%%%%%%%%%%
    objf  = [];
    objfn = {};
    ndata = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;
            

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues1 = simulate_DR_IPTG_TetR(para,data,ParaNames,model);
 
            
        %plot

            figure(1)
            %subplot(2,3,1);
            hold on;
            plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2);
            
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('LacI (1st repression coefficient) & TetRTup1', 'FontSize',20)
            %legend({'simulation' 'data'}, 'Location', 'southeast')


        %% module 2: PAct1_LacI(W220F)_tCyc1
            para = paraValues;
            

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep2', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues2 = simulate_DR_IPTG_TetR(para,data,ParaNames,model);
            
        %plot

            figure(2)
            %subplot(2,3,2);
            hold on;
            plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2);
            %errorbar(log10(data.dose),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence', 'FontSize', 18)
            title('LacI (2nd repression coefficient) & TetRTup1', 'FontSize', 20)
            %legend({'simulation-W220F' 'data-W220F'}, 'Location', 'southeast')


        %% module 3: P4Lacn.2_LacI(W220F,Q60G,T167A)_tCyc1
            para = paraValues;
      
            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep3', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues3 = simulate_DR_IPTG_TetR(para,data,ParaNames,model);
            
        %plot

            figure(3)
            %subplot(2,3,3);
            hold on;
            plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2);
            %errorbar(log10(data_W220F_Q60G_T167A.dose),DataMeans,DataStd,'o');
            xlabel('log IPTG (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence', 'FontSize', 18)
            title('LacI (3rd repression coefficient) & TetRTup1', 'FontSize', 20)
            %legend({'simulation-W220F-Q60G-T167A' 'data-W220F-Q60G-T167A'}, 'Location', 'southeast')

figure(4)
hold on;
plot(log10(data.dose),SimFluoValues1,'-', 'LineWidth', 2, 'Color','red');
plot(log10(data.dose),SimFluoValues2,'-', 'LineWidth', 2, 'Color', 'blue');
plot(log10(data.dose),SimFluoValues3,'-', 'LineWidth', 2, 'Color', 'green');
xlabel('log IPTG (nM)', 'FontSize', 18)
ylabel('mean Fluorescence', 'FontSize', 18)
title('LacI (all repression coefficient) & TetRTup1', 'FontSize', 20)
hold off
legend('LacI (1st repression coefficient) & TetRTup1','LacI (2nd repression coefficient) & TetRTup1','LacI (3rd repression coefficient) & TetRTup1', 'Location', 'northeastoutside')



end

        

