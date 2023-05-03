function SimFluoValues = simulate_DR_IPTG(para, data_IPTG,data, ParaNames,model)
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
    SimFluoValuesAll = cell(height(dose), 1);
    parfor kdose = 1:height(dose)
        d1_copy = d1; % create a copy of d1 for this iteration
        d1_copy.Amount = data_IPTG.dose(kdose,1); % modify the copy
        
        SimFluoValues = zeros(height(dose_aTc), 1);
        
        parfor adose = 1:height(dose_aTc)
            d2_copy = d2; % create a copy of d2 for this iteration
            d2_copy.Amount = data.dose(adose,1); % modify the copy
            
            try      
                [t,sd,species] = sbiosimulate(model.mw_sbmod1, [d1_copy, d2_copy]);
                index = ismember(species, 'Citrine');
                SimCitrineValues = sd(end,index);
                SimFluoValues(adose) = x_scal_aTc*SimCitrineValues + data.empty;
            catch ME
                disp(ME)
            end
        end
        
        SimFluoValuesAll{kdose} = SimFluoValues;
    end
    
    % Concatenate results from all iterations
    SimFluoValues = cat(2, SimFluoValuesAll{:});
end