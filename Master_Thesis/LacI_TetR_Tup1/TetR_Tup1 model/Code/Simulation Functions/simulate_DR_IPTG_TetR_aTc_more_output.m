function [SimFluoValues, time_IPTG, time_aTc, SimFluoValues_over_time_aTc, SimFluoValues_over_time_IPTG] = simulate_DR_IPTG_TetR_aTc_more_output(para, data_IPTG,data, ParaNames,model)
% simulates the fluorescent values with parameters para using the
% timepoints from the data


% Set simulation Time
    tmax = data_IPTG.time;
    tmax_aTc = data.time;
    %index = ismember(ParaNames, 'indTime');
    %indTime = para(index);
    indTime = 2000; % Alternative line as indTime is not a parameter in the test model
    cs = getconfigset(model.mw_sbmod1);
    cs.StopTime = (tmax+indTime+tmax_aTc);
    cs.SolverType = 'sundials';
    
% Doses, Scaling Factor
    dose = data_IPTG.dose;
    d1 = sbiodose('d1','schedule');
    d1.Amount = 0;
    d1.time = indTime;
    d1.TargetName = 'IPTG';
    d1.Active = true;
    %x_scal = 0.0077*(data.tdh3 - data.empty); % scaling factor
    x_scal =  (data_IPTG.tdh3 - data_IPTG.empty)*0.0077*(0.0173+0.0077)/0.0173;  % scaling factor (Fact - F0)*rho*(km+rho)/km       

    dose_aTc = data.dose;
    d2 = sbiodose('d2', 'schedule');
    d2.Amount = 0;
    d2.time = indTime;
    d2.TargetName = 'aTc';
    d2.Active = true;
    x_scal_aTc = (data.tdh3 - data.empty)*0.0077*(0.0173+0.0077)/0.0173;

% Apply parametervalues to the model
    for index = 1:length(para)
        model.mw_sbmod1.Parameters(index).Value = para(index,1);
    end
    
% Simuate model
    time_aTc = cell(height(dose_aTc), 1);
    time_IPTG = cell(height(dose),1);
    SimFluoValues_over_time_aTc = cell(height(dose_aTc),1);
    SimFluoValues_over_time_IPTG = cell(height(dose),1);
    SimFluoValues = zeros(height(dose),1);
    for kdose = 1:height(dose)
        d1.Amount = data_IPTG.dose(kdose,1); % IPTG To Do: remove nMperUnit, place it inside the model
        for adose = 1:height(dose_aTc)
            d2.Amount = data.dose(adose,1);
        try      
            [t,sd,species] = sbiosimulate(model.mw_sbmod1, [d1, d2]);
            index = ismember(species, 'Citrine');
            time_aTc{adose} = t;
            SimCitrineValues = sd(end,index);
            SimFluoValues(adose,kdose) = x_scal_aTc*SimCitrineValues + data.empty;
            SimValues_over_time_aTc = sd(:,index);
            SimFluoValues_over_time_aTc{adose} = x_scal_aTc.*SimValues_over_time_aTc + data.empty;

        catch ME
            disp(ME)
        end
        
        end 

        time_IPTG{kdose} = t;
        SimValues_over_time_IPTG = sd(:,index);
        SimFluoValues_over_time_IPTG{kdose} = x_scal_aTc.*SimValues_over_time_IPTG + data.empty;

    end  