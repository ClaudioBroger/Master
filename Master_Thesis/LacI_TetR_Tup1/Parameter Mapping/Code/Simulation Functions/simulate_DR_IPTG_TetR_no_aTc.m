function [SimFluoValues, SimFluoValues_time, time_IPTG] = simulate_DR_IPTG_TetR_no_aTc(para, data, ParaNames,model)
% simulates the fluorescent values with parameters para using the
% timepoints from the data


% Set simulation Time
    tmax = data.time;
    indTime = 2000; % Alternative line as indTime is not a parameter in the test model
    cs = getconfigset(model.mw_sbmod1);
    cs.StopTime = (tmax+indTime);
    cs.SolverType = 'sundials';
    
% Doses, Scaling Factor
    dose = data.dose;
    d1 = sbiodose('d1','schedule');
    d1.Amount = 0;
    d1.time = indTime;
    d1.TargetName = 'IPTG';
    d1.Active = true;

    x_scal =  (data.tdh3 - data.empty)*0.0077*(0.0173+0.0077)/0.0173;  % scaling factor (Fact - F0)*rho*(km+rho)/km       

% Apply parametervalues to the model
    for index = 1:length(para)
        model.mw_sbmod1.Parameters(index).Value = para(index,1);
    end
    
% Simuate model
    time_IPTG = cell(height(dose),1);
    SimFluoValues_time = cell(height(dose),1);
    SimFluoValues = zeros(height(dose),1);
    for kdose = 1:height(dose)
        d1.Amount = data.dose(kdose,1); % IPTG To Do: remove nMperUnit, place it inside the model  
        try      
            [t,sd,species] = sbiosimulate(model.mw_sbmod1, d1);
            index = ismember(species, 'Citrine');
            SimCitrineValues = sd(end,index);
            SimFluoValues(kdose) = x_scal*SimCitrineValues + data.empty;
            SimFluoValues_over_time = sd(:,index)';

        catch ME
            disp(ME)
        end
        SimFluoValues_time{kdose} = x_scal.*SimFluoValues_over_time + data.empty;
        time_IPTG{kdose} = t;
    end     
end  