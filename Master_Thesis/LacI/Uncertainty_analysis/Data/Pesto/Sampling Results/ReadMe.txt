These results are obtained using the PESTO toolbox.


some usefull code: 
confLevels = [0.95];
parameters = getParameterConfidenceIntervals(parameters, confLevels, optionsPesto);

samplingPlottingOpt = PestoPlottingOptions();
samplingPlottingOpt.S.plot_type = 1;
samplingPlottingOpt.S.ind = 1;
h = figure('Name','plotParameterSamples - 1D');
plotParameterSamples(parameters,'1D',fh,[],samplingPlottingOpt);
