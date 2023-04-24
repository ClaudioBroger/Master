%PESTO uncertainty

load('1704_fitModel_Pesto.mat')

dataPath_aTc = ('/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 model/');
%load model for uncertainty analysis
model = sbioloadproject('TetR_Tup1_model');
FlPlot = true;
num_draws = 30;
counter = 0;

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));


parameters.name = cell2table(parameters.name);

[~,index] = ismember(ParaNames,parameters.name{:,:});
index = nonzeros(index);

min_par = parameters.min(index);
max_par = parameters.max(index);


prob = parameters.S.par(index,:)';

for draw=1:num_draws
    for k=1:length(min_par)
        pd = fitdist(prob(:,k),'Kernel');
         X = min_par(k):0.001:max_par(k);
         Y = pdf(pd,X);
%          figure(1)
%          subplot(2,6,k)
         plot(Y)
         rand_parameter(draw,k) = randsample(X,1,true,Y);
         
    end
end



% num_draws = 30;
% for k = 1:num_draws
%     for j = 1:length(CI)
%         rand_parameter(k,j) = CI(j,1) + (CI(j,2) - CI(j,1)).*rand;
%     end
% end



rand_parameter = array2table(rand_parameter);

rand_parameter.Properties.VariableNames(1:5) = ParaNames(1:5);
rand_parameter.Properties.VariableNames(6:13) = ParaNames(9:16);

rand_parameter.mu(:) = 0.0077;
rand_parameter.kmaturation(:) = 0.0173;


rand_parameter.nMperUnit(:) = 10^1.1;
rand_parameter.a(:) = 1;
rand_parameter.f(:) = 1;
rand_parameter.g(:) = 1;
rand_parameter.b(:) = 1;
[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));
rand_parameter = rand_parameter(:,[string(ParaNames)]);
rand_parameter.dCit(:) = 0;

rand_parameter(:,1:3) = array2table(10.^table2array(rand_parameter(:,1:3)));
rand_parameter(:,5) = array2table(10.^table2array(rand_parameter(:,5)));
rand_parameter(:,10:15) = array2table(10.^table2array(rand_parameter(:,10:15)));
rand_parameter.p1(:) = 5;
rand_parameter.p2(:) = 5;
save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI/Uncertainty_analysis/Data/Pesto/',datestr(now, 'dd-mmm-yyyy'),'_rand_parameter_Pesto.mat'],'rand_parameter')

C = linspecer(num_draws);
for draw = 1:num_draws
    
    %define parameter values as the values randomly drawn
    paraValues = table2array(rand_parameter(draw,:));
    paraValues = paraValues'; 
    %TFlMode = FlMode;
    
     
        % INITIALIZE OUTPUTS %
        %%%%%%%%%%%%%%%%%%%%%%
        objf  = [];
        objfn = {};
        ndata = 0;
    
        %% module 1: ptetTetR + ptetCitrine
            para = paraValues;

            dataPos_aTc = strcat(dataPath_aTc, "dataSRtup1.mat");
            data_aTc = load(dataPos_aTc);
            
            
            
            NamestoZero = setdiff(ParaNames,{'kactTetR','k7tetCit', 'kL7tet', 'dTetR', 'dCit', 'thetaTetR', 'nTetR', 'f', 'nTup1', 'KdTetR', 'growthFix', 'nMperUnit', 'kmaturation', 'indTime' , 'thetaTup1', 'a', 'p1', 'p2'});                            
            IdxToZero = find(ismember(ParaNames, NamestoZero)) ;           
            para(IdxToZero) = 0;

            SimFluoValues1 = simulate_DR_aTc(para,data_aTc,ParaNames,model);
            DataMeans = data_aTc.data.means;
            DataStd = data_aTc.data.std;
            data_aTc.dose(1) = 1;
            

%         %plot
        if FlPlot
            figure(2)
            hold on;
            plot(log10(data_aTc.dose),SimFluoValues1,'-', 'LineWidth', 2, 'DisplayName',strcat('Simulation-rep1-', num2str(draw)), 'Color', C(draw,:));
            %errorbar(log10(data.dose),DataMeans,DataStd,'o', 'HandleVisibility','off');
            xlabel('log aTc (nM)', 'FontSize', 18)
            ylabel('mean Fluorescence','FontSize', 18)
            title('P4Lacn.2-cit + PAct1-LacI, 1st repression coefficient', 'FontSize',20)
            
        end
        hold off
        legend("show", 'Location', 'northeastoutside')
%     
            % From here on the steps are repeated as in the module (with
            % first repression coefficient) above
        

end
