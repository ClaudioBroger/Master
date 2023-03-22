%Load data
load('2023_02_09_fitModel_Hyperspace_npoints20_Eline.mat')
%Load model
model = sbioloadproject('LacImodel');
%Load model settings
ModelSettings_4;
%Specify number of random draws
num_draws = 200;

%find estimated parameters in data
paramSpecs = paramSpecs(Settings.model.PIdx,:);
if ~isfield(OutV,'colnames')
    OutV.colnames = paramSpecs.names;
end
NotProjectedIdx = find(ismember(paramSpecs.names,OutV.colnames));
% paraoptestimated = paraopt(NotProjectedIdx); % to add fitted parameter point

OutV.V = array2table(OutV.V, 'VariableNames', OutV.colnames);
OutV.V(:,7) = {1.4};

%create table
paramSpecsestimated = paramSpecs(NotProjectedIdx,:);
data = OutV;
if isfield(OutV,'V')
    data.rowmat = data.V;
end

%take minimum and maximum from estimated parameters
paramNames = paramSpecsestimated.names;
bmin = paramSpecsestimated.bmin;
bmax = paramSpecsestimated.bmax;
bmin(paramSpecsestimated.islog==1) = log10(bmin(paramSpecsestimated.islog==1));
bmax(paramSpecsestimated.islog==1) = log10(bmax(paramSpecsestimated.islog==1));
bmin(7) = 1.4;
bmax(7) = 1.4;
paramtoplot = paramNames;
n = length(paramtoplot);

%%%nMperUnit was at position 7 -> need to take out column 7 from
%%%data.rowmat to not include it in further steps
%data.rowmat(:,7) = [];

%create minimum and maximum vectors from data
%minimum and maximum vectors not the same length -> fill the rest with 0s
indices = zeros(1,n);
x = zeros(height(data.rowmat(:,1)),n);
xmin = zeros(1,n);
xmax = zeros(1,n);
x = array2table(x);
for k =1:n
    indices(k) = find(ismember(paramNames,paramtoplot{k}));
    x(:,k) =  data.rowmat(:,k);
    xmin(k) = bmin(indices(k));
    xmax(k) = bmax(indices(k));
end
x = table2array(x);

%Fit the probability density function depending on the parameter values
%draw random values from the probability density function
%many 0s -> onyl chose random parameters which are not 0
for draw=1:num_draws
    for k=1:n
        pd = fitdist(x(:,k),'Kernel');
         X = xmin(k):0.001:xmax(k);
         Y = pdf(pd,X);
         rand_parameter(draw,k) = randsample(X,1,true,Y);
         
    end
end
rand_parameter = array2table(rand_parameter, 'VariableNames',paramtoplot);

%save random parameters
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/',datestr(now, 'dd-mmm-yyyy'),'_rand_parameter_ec50.mat'],'rand_parameter')



