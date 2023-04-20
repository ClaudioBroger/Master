
load('0324_PestoResults_cleanup.mat')

for k = 1:17
    figure(1)
    subplot(3,6,k)
    histogram(parameters.S.par(k,:))
end

load('2023_02_09_fitModel_Hyperspace_npoints20_Eline.mat')

for k = 1:17
    figure(2)
    subplot(3,6,k)
    histogram(OutV.V(:,k))
end

%%Parameter distribution PESTO before drawing
load('0324_PestoResults_cleanup.mat')

names = ["PAct1LacI", "P4LacnCit", "dLacI", "LacIrepWT", "KdLacI", "nLacI", "nMperUnit", "LacIrepW220F", "P4LacnLacI", "P4LacnLacIL", "LacIrep3mut", "pt7LacI", "P3Lacn5Cit", "P3Lacn5CitL", "dLacIpt7", "nLacIP3", "LacIrep3mutP3"]';

parameters.name(1) = {"PAct1LacI"};
parameters.name(2) = {"P4LacnLacI"};
parameters.name(3) = {"pt7LacI"};
parameters.name(4) = {"P3Lacn5Cit"};
parameters.name(5) = {"P4LacnCit"};
parameters.name(6) = {"dLacI"};
parameters.name(7) = {"dLacIpt7"};
parameters.name(8) = {"P4LacnLacIL"};
parameters.name(9) = {"P3Lacn5CitL"};
parameters.name(10) = {"LacIrepWT"};
parameters.name(11) = {"LacIrepW220F"};
parameters.name(12) = {"LacIrep3mut"};
parameters.name(13) = {"LacIrep3mutP3"};
parameters.name(14) = {"KdLacI"};
parameters.name(15) = {"nLacI"};
parameters.name(16) = {"nLacIP3"};
parameters.name(17) = {"nMperUnit"};

dataPath = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/DataForEstimation/');
%load model for uncertainty analysis
model = sbioloadproject('LacImodel_uncertainty');
FlPlot = true;
num_draws = 30;
counter = 0;


for k = 1:length(names)
    index(k) = find(ismember(string(parameters.name), names(k)));
end

parameters.min = parameters.min(index);
parameters.max = parameters.max(index);

min_par = parameters.min;
max_par = parameters.max;

CI = parameters.CI.S(index,:);

prob = parameters.S.par(index,:)';


    for k=1:length(CI)
        pd = fitdist(prob(:,k),'Kernel');
         X = min_par(k):0.001:max_par(k);
         Y = pdf(pd,X);
         figure(1)
         subplot(3,6,k)
         plot(X,Y,'Color','black','LineWidth',2);
         axis 'auto xy'
         title(names(k), 'FontSize', 12);
         ax = gca;
         ax.FontSize = 20;
         
    end

%%Parameterdistribution hyperspace after drawing
load("30-Mar-2023_rand_parameter_200.mat")
rand_parameter.nMperUnit = [];
for k = 1:width(rand_parameter)
    subplot(3,6,k)
    histogram(table2array(rand_parameter(:,k)));
    axis 'auto xy'
    title(names(k), 'FontSize', 12);
    ax = gca;
    ax.FontSize = 20;
end

%%Parameterdistribution PESTO after drawing
load('03-Apr-2023_rand_parameter_Pesto.mat')
rand_parameter.dCit = [];
rand_parameter.mu = [];
rand_parameter.kmaturation = [];
rand_parameter.nMperUnit = [];
for k = 1:width(rand_parameter)
    subplot(3,3,k)
    histogram(table2array(rand_parameter(:,k)));
    axis 'auto xy'
    title(rand_parameter.Properties.VariableNames(k), 'FontSize', 12);
    ax = gca;
    ax.FontSize = 20;
end

%%Parameterdistribution PESTO CI after drawing
load('03-Apr-2023_rand_parameter_Pesto_cI-200.mat')
rand_parameter.dCit = [];
rand_parameter.mu = [];
rand_parameter.kmaturation = [];
rand_parameter.nMperUnit = [];
for k = 1:width(rand_parameter)
    subplot(3,3,k)
    histogram(table2array(rand_parameter(:,k)));
    axis 'auto xy'
    title(rand_parameter.Properties.VariableNames(k), 'FontSize', 12);
    ax = gca;
    ax.FontSize = 20;
end
