function [growth_rate_IPTG, growth_rate_aTc, time_aTc, growth_rate, SimTetRTup1Values_aTc, SimTetRTup1freeValues_aTc, SimLacIValues_aTc, SimLacIfreeValues_aTc, SimTetRTup1Values_IPTG, SimTetRTup1freeValues_IPTG, SimLacIValues_IPTG, SimLacIfreeValues_IPTG, time_IPTG] = simulate_DR_IPTG_TetR_aTc_PMA1(para, data_IPTG,data, ParaNames,model)
% simulates the fluorescent values with parameters para using the
% timepoints from the data


% Set simulation Time
    tmax = data_IPTG.time;
    %index = ismember(ParaNames, 'indTime');
    %indTime = para(index);
    indTime = 2000; % Alternative line as indTime is not a parameter in the test model
    cs = getconfigset(model.mw_sbmod1);
    cs.StopTime = (tmax+indTime);
    cs.SolverType = 'sundials';
     set(cs.SolverOptions, 'MaxStep', 10);
% Initialize configset for analysis run.
    

% Doses, Scaling Factor
    dose = data_IPTG.dose;
    d1 = sbiodose('d1','schedule');
    d1.Amount = 0;
    d1.time = indTime;
    d1.TargetName = 'IPTG';
    d1.Active = true;
    
    dose_aTc = data.dose;
    d2 = sbiodose('d2', 'schedule');
    d2.Amount = 0;
    d2.time = indTime;
    d2.TargetName = 'aTc';
    d2.Active = true;
    
% Apply parametervalues to the model
    for index_PMA1 = 1:length(para)
        model.mw_sbmod1.Parameters(index_PMA1).Value = para(index_PMA1,1);
    end
    
% Simuate model
     growth_rate_IPTG = cell(height(dose), 1);
     SimPMA1Values = zeros(height(dose),1);
     growth_rate_aTc = cell(height(dose_aTc), 1);
     time_aTc = cell(height(dose_aTc), 1);
     time_IPTG = cell(height(dose),1);
     SimTetRTup1Values_aTc = cell(height(dose_aTc),1);
     SimTetRTup1freeValues_aTc = cell(height(dose_aTc),1);
     SimLacIValues_aTc = cell(height(dose_aTc),1);
     SimLacIfreeValues_aTc = cell(height(dose_aTc),1);
     SimTetRTup1Values_IPTG = cell(height(dose),1);
     SimTetRTup1freeValues_IPTG = cell(height(dose),1);
     SimLacIValues_IPTG = cell(height(dose),1);
     SimLacIfreeValues_IPTG = cell(height(dose),1);
    for kdose = 1:height(dose)
        d1.Amount = data_IPTG.dose(kdose,1); % IPTG To Do: remove nMperUnit, place it inside the model
        for adose = 1:height(dose_aTc)
            d2.Amount = data.dose(adose,1);
        try      
            [t,sd,species] = sbiosimulate(model.mw_sbmod1, [d1, d2]);
            index_PMA1 = ismember(species, 'PMA1');
            SimPMA1Values = sd(end,index_PMA1);
            SimPMA1Values_over_time = sd(:,index_PMA1)';
            growth_rate(adose,kdose) = max(0,para(22)  + (para(21)  - para(22) )/(1 + exp(para(20)*(log((SimPMA1Values)*para(14))-para(19)))));
             growth_rate_aTc{adose} = max(0, para(22)  + (para(21)  - para(22) )./(1 + exp(para(20)*(log((SimPMA1Values_over_time)*para(14))-para(19)))));
            time_aTc{adose} = t;
            index_TetRTup1 = ismember(species, 'TetRTup1');
            SimTetRTup1Values_aTc{adose} = sd(:, index_TetRTup1);
            index_TetRTup1free = ismember(species, 'TetRTup1free');
            SimTetRTup1freeValues_aTc{adose} = sd(:, index_TetRTup1free);
            index_LacI = ismember(species, 'LacI');
            SimLacIValues_aTc{adose} = sd(:, index_LacI);
            index_LacIfree = ismember(species, 'LacIfree');
            SimLacIfreeValues_aTc{adose} = sd(:, index_LacIfree);

            
        catch ME
            disp(ME)
        end
        
        end

        growth_rate_IPTG{kdose} = max(0,para(22)  + (para(21)  - para(22) )./(1 + exp(para(20)*(log((SimPMA1Values_over_time)*para(14))-para(19)))));
        time_IPTG{kdose} = t;
        SimTetRTup1Values_IPTG{kdose} = sd(:,index_TetRTup1);
        SimTetRTup1freeValues_IPTG{kdose} = sd(:,index_TetRTup1free);
        SimLacIValues_IPTG{kdose} = sd(:,index_LacI);
        SimLacIfreeValues_IPTG{kdose} = sd(:,index_LacIfree);


    end  