function SimFluoValues = simulateDRaTc(para, data, ParaNames,model)
% simulates the fluorescent values with parameters para using the
% timepoints from the data
    tmax = data.time;
    index = ismember(ParaNames, 'indTime');
    indTime = para(index);
    tspan = [0:1:(tmax+indTime)];   
    dose = data.dose;
    %x_scal = 0.0077*(data.tdh3 - data.empty); % scaling factor
    x_scal =  (data.tdh3 - data.empty)*0.0077*(0.0173+0.0077)/0.0173;  % scaling factor (Fact - F0)*rho*(km+rho)/km       
    fnhandle = model;


    SimFluoValues = zeros(length(dose),1);
    for kdose = 1:length(dose)
        indexaTc = ismember(ParaNames, 'atcAdded');
        para(indexaTc) = dose(kdose); %aTcAdded   
        try
            output =  fnhandle(tspan, [], para);  
            index = ismember(output.states, 'Citrine');
            SimCitrineValues = output.statevalues(end,index);
            SimFluoValues(kdose) = x_scal*SimCitrineValues + data.empty;
           
        catch ME
            disp(ME)
        end
    end     
end  
