load('1704_fitModel_Pesto.mat')

num_draws = 30;

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));

parameters.name(7) = {"kLacI"};
parameters.name(3) = {"kTetRTup1"};
parameters.name(13) = {"dTetRTup1"};
parameters.name(14) = {"degtag"};
parameters.name(23) = {"LacIrep"};
parameters.name(25) = {"LacIrep2"};
parameters.name(26) = {"LacIrep3"};
parameters.name(5) = {"kCit"};
parameters.name(17) = {"dCit"};
parameters.name(18) = {"CitL"};
parameters.name(20) = {"TetRTup1L"};
parameters.name(31) = {"nTetRTup1"};
parameters.name(22) = {"TetRTup1rep"};
parameters.name(2) = {"kTetR"};
parameters.name(19) = {"LacIL"};


parameters.name = cell2table(parameters.name);

[~,index] = ismember(ParaNames,parameters.name{:,:});
index = nonzeros(index);

min_par = parameters.min(index);
max_par = parameters.max(index);

prob = parameters.S.par(index,:)';

for draw=1:num_draws
    for k=1:width(prob)
        pd = fitdist(prob(:,k),'Kernel');
         X = min_par(k):0.001:max_par(k);
         Y = pdf(pd,X);
%          figure(1)
%          subplot(2,6,k)
         plot(Y)
         rand_parameter(draw,k) = randsample(X,1,true,Y);
         
    end
end

rand_parameter(:,1:12) = 10.^rand_parameter(:,1:12);
rand_parameter(:,16) = 10.^rand_parameter(:,16);
rand_parameter(:,18) = 10.^rand_parameter(:,18);
rand_parameter(:,19:21) = 10.^rand_parameter(:,19:21);

rand_parameter = array2table(rand_parameter);
rand_parameter.Properties.VariableNames(1:16) = ParaNames(1:16);
rand_parameter.Properties.VariableNames(17) = ParaNames(18);
rand_parameter.Properties.VariableNames(18) = ParaNames(20);
rand_parameter.Properties.VariableNames(19:21) = ParaNames(21:23);
rand_parameter.nMperUnit(:) = 10^1.1;

rand_parameter.mu(:) = 0.0077;
rand_parameter.kmaturation(:) = 0.0173;
rand_parameter.dCit(:) = 0;

rand_parameter = rand_parameter(:,[string(ParaNames)]);

save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/TetR_Tup1 NF/Results/',datestr(now, 'dd-mmm-yyyy'),'_rand_parameter_LacI_NF_TetR.mat'],'rand_parameter')

