C = linspecer(height(data.dose));
counter = 0;
for draw = 1:height(rand_parameter)
    
    %define parameter values as the values randomly drawn
     paraValues = table2array(rand_parameter(draw,:));
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
            SimFluoValues1(height(data.dose)+1:end,:) = [];



            objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep2', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep', 'kTetR', 'dTetR'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues2 = simulate_DR_IPTG_TetR_aTc(para,data_IPTG,data,ParaNames,model);
            SimFluoValues2(height(data.dose)+1:end,:) = [];


            objf  = [];
        objfn = {};
        ndata = 0;

        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            NamestoZero = setdiff(ParaNames,{'kLacI','dLacI', 'kTetRTup1', 'dTetRTup1', 'degtag','LacIrep3', 'kCit', 'CitL', 'TetRTup1L', 'nLacI', 'indTime', 'nTetRTup1', 'KdLacI', 'KdTetR', 'kmaturation', 'nMperUnit', 'mu', 'TetRTup1rep', 'kTetR', 'dTetR'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;


            SimFluoValues3 = simulate_DR_IPTG_TetR_aTc(para,data_IPTG,data,ParaNames,model);
            SimFluoValues3(height(data.dose)+1:end,:) = [];



%% Combine SimFluoValues files from all draws

SimFluoValues1 = array2table(SimFluoValues1);
SimFluoValues2 = array2table(SimFluoValues2);
SimFluoValues3 = array2table(SimFluoValues3);

SimFluoValues1.draw(:) = draw;
SimFluoValues2.draw(:) = draw;
SimFluoValues3.draw(:) = draw;

SimFluoValues1.rep(:) = rand_parameter.LacIrep(draw);
SimFluoValues2.rep(:) = rand_parameter.LacIrep2(draw);
SimFluoValues3.rep(:) = rand_parameter.LacIrep3(draw);

SimFluoValues1.aTc_dose = data.dose;
SimFluoValues2.aTc_dose = data.dose;
SimFluoValues3.aTc_dose = data.dose;

if draw == 1
    SimFluoValues1_combined_aTc = SimFluoValues1;
    SimFluoValues2_combined_aTc = SimFluoValues2;
    SimFluoValues3_combined_aTc = SimFluoValues3;
else

    SimFluoValues1_combined_aTc = [SimFluoValues1_combined_aTc; SimFluoValues1];
    SimFluoValues2_combined_aTc = [SimFluoValues2_combined_aTc; SimFluoValues2];
    SimFluoValues3_combined_aTc = [SimFluoValues3_combined_aTc; SimFluoValues3];
end
end

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/Results/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues1_combined_LacI_TetRTup1', '.mat'], 'SimFluoValues1_combined_aTc');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/Results/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues2_combined_LacI_TetRTup1', '.mat'], 'SimFluoValues2_combined_aTc');
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/Results/', datestr(now, 'dd-mmm-yyyy'),'SimFluoValues3_combined_LacI_TetRTup1', '.mat'], 'SimFluoValues3_combined_aTc');







