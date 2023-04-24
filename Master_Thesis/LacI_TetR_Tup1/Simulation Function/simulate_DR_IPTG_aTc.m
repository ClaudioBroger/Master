function SimFluoValues = simulate_DR_IPTG(para, data_IPTG,data_aTc, ParaNames,model)
% simulates the fluorescent values with parameters para using the
% timepoints from the data


% Set simulation Time
    tmax = data_IPTG.data.time;
    %index = ismember(ParaNames, 'indTime');
    %indTime = para(index);
    indTime = 2000; % Alternative line as indTime is not a parameter in the test model
    cs = getconfigset(model.mw_sbmod1);
    cs.StopTime = (tmax+indTime);
    cs.SolverType = 'sundials';
    
% Doses, Scaling Factor
    dose_IPTG = data_IPTG.dose;
    d1 = sbiodose('d1','schedule');
    d1.Amount = 0;
    d1.time = indTime;
    d1.TargetName = 'IPTG';
    d1.Active = true;

    dose_aTc = data_aTc.data.dose;
    d2 = sbiodose('d2','schedule');
    d2.Amount = 0;
    d2.time = indTime;
    d2.TargetName = 'aTc';
    d2.Active = true;
    %x_scal = 0.0077*(data.tdh3 - data.empty); % scaling factor
    x_scal =  (data_aTc.data.tdh3 - data_aTc.data.empty)*0.0077*(0.0173+0.0077)/0.0173;  % scaling factor (Fact - F0)*rho*(km+rho)/km       

% Apply parametervalues to the model
    for index = 1:length(para)
        model.mw_sbmod1.Parameters(index).Value = para(index,1);
    end
    
% Simuate model

    SimFluoValues = zeros(height(dose_IPTG),1);
    for kdose = 1:height(dose_IPTG)
        d1.Amount = data_IPTG.dose(kdose,1);% IPTG To Do: remove nMperUnit, place it inside the model  
    for adose = 1:height(dose_aTc)
        d2.Amount = data_aTc.data.dose(adose,1);
        try      
            [t,sd,species] = sbiosimulate(model.mw_sbmod1, [d1, d2]);
            index = ismember(species, 'Citrine');
            SimCitrineValues = sd(end,index);
            SimFluoValues(kdose) = x_scal*SimCitrineValues + data.data.empty;
        catch ME
            disp(ME)
        end
        
    end     
end  