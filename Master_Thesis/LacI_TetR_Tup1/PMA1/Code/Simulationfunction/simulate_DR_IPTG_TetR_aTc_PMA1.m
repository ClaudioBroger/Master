function growth_rate = simulate_DR_IPTG_TetR_aTc(para, data_IPTG,data, ParaNames,model)
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

    SimPMA1Values = zeros(height(dose),1);
    growth_rate_o_t = zeros(height(dose),1);
    for kdose = 1:height(dose)
        d1.Amount = data_IPTG.dose(kdose,1); % IPTG To Do: remove nMperUnit, place it inside the model
        for adose = 1:height(dose_aTc)
            d2.Amount = data.dose(adose,1);
        try      
            [t,sd,species] = sbiosimulate(model.mw_sbmod1, [d1, d2]);
            index = ismember(species, 'PMA1');
            SimPMA1Values = sd(end,index);
            SimPMA1Values_over_time = sd(:,index)';
            %SimPMA1Values(adose, kdose) = x_scal_aTc*SimPMA1Values + data.empty;
            growth_rate(adose,kdose) = max(0,para(29)  + (para(28)  - para(29) )/(1 + exp(para(27)*(log((SimPMA1Values)*para(17))-para(26)))));
%             growth_rate_o_t(adose,:) = max(0,para(29)  + (para(28)  - para(29) )./(1 + exp(5*(log((SimPMA1Values_over_time)*para(17))-5))));
%             time = t;
        catch ME
            disp(ME)
        end
        
        end   
%     if kdose == 1    
%         growth_rate_over_time = growth_rate_o_t;
%     else
%         growth_rate_over_time = [growth_rate_o_t; growth_rate_over_time];
%     end
    end  