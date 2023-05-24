C = linspecer(height(data.dose));
counter = 0;
% rand_parameter.kTetR = [];
% rand_parameter.dTetR = [];


for k = 1:width(rand_parameter)
    special2 = rand_parameter(2,:);
    special2 = table2array(special2);
    special2(k) = special2(k)*0,2;
    paraValues = special2;
     paraValues = paraValues';
    %TFlMode = FlMode;
    
     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;
    
        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep', 'kTetR', 'dTetR'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

     

            SimFluoValues1 = simulate_DR_IPTG_TetR_aTc(para,data_IPTG,data,ParaNames,model);
           
            figure(k)
            for a = 1:height(data.dose)
                subplot(3,1,1)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(data_IPTG.dose),SimFluoValues1(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence','FontSize', 18)
                title('LacI (1st repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end
            
            legend("show", 'Location', 'northeastoutside')

%         %plot

            
            
        
        
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

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep2', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep', 'kTetR', 'dTetR'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues2 = simulate_DR_IPTG_TetR_aTc(para,data_IPTG,data,ParaNames,model);
            figure(k)
            for a = 1:height(data.dose)
                subplot(3,1,2)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(data_IPTG.dose),SimFluoValues2(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence','FontSize', 18)
                title('LacI (2nd repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end
            
            legend("show", 'Location', 'northeastoutside')
           
            
        %plot

            

     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep3', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep', 'kTetR', 'dTetR'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues3 = simulate_DR_IPTG_TetR_aTc(para,data_IPTG,data,ParaNames,model);



        %plot

            figure(k)
            for a = 1:height(data.dose)
                subplot(3,1,3)
                hold on;
                
                % Round the data.dose(a) value to 2 decimal places
                rounded_dose = round(data.dose(a) * 100) / 100;
                
                % Convert the rounded value to a string and concatenate with the label
                plot(log10(data_IPTG.dose),SimFluoValues3(a,:),'-', 'LineWidth', 2, 'DisplayName', strcat('aTc dose: ', num2str(rounded_dose), ' nM'), 'Color', C(a,:));
                
                xlabel('log IPTG (nM)', 'FontSize', 18)
                ylabel('mean Fluorescence','FontSize', 18)
                title('LacI (3rd repression coefficient) & TetRTup1', 'FontSize',20)
                
            
            end
            
            legend("show", 'Location', 'northeastoutside')
            objf  = [];
            objfn = {};
            ndata = 0;
end