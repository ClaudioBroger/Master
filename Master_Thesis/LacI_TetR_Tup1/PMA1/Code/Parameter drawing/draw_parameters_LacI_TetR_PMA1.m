load('1704_fitModel_Pesto.mat')

num_draws = 30;

[ParaNames] = (get(model.mw_sbmod1.parameters, {'Name'}));

parameters.name(6) = {"kLacI"};
parameters.name(3) = {"kTetRTup1"};
parameters.name(13) = {"dTetRTup1"};
parameters.name(14) = {"degtag"};
parameters.name(23) = {"LacIrep"};
parameters.name(25) = {"LacIrep2"};
parameters.name(26) = {"LacIrep3"};
parameters.name(5) = {"kCit"};
parameters.name(17) = {"dCit"};
parameters.name(18) = {"PMA1L"};
parameters.name(20) = {"TetRTup1L"};
parameters.name(31) = {"nTetRTup1"};
parameters.name(22) = {"TetRTup1rep"};
parameters.name(2) = {"kTetR"};
parameters.name(4) = {"kPMA1"};
parameters.name(21) = {"TetRrep"};
parameters.name(11) = {"kLacTetRTup1"}




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

rand_parameter(:,1:11) = 10.^rand_parameter(:,1:11);
rand_parameter(:,15) = 10.^rand_parameter(:,15);
rand_parameter(:,17:22) = 10.^rand_parameter(:,17:22);
rand_parameter(:,24) = 10.^rand_parameter(:,24);
rand_parameter(:,26) = 10.^rand_parameter(:,26);
rand_parameter = array2table(rand_parameter);
rand_parameter.Properties.VariableNames(:) = table2array(parameters.name(index,:));

rand_parameter.nMperUnit(:) = 10^1.1;

rand_parameter.mu(:) = 0.0077;
rand_parameter.kmaturation(:) = 0.0173;
rand_parameter.dCit(:) = 0;
rand_parameter.p1(:) = 5;
rand_parameter.p2(:) = 5;
rand_parameter.growthMIN(:) = 0;
rand_parameter.growthMAX(:) = 0.0077;
rand_parameter.f(:) = 1;
rand_parameter.g(:) = 1;
rand_parameter.PMAfactor(:) = 1;

rand_parameter = rand_parameter(:,[string(ParaNames)]);


save(['/Users/claudiobroger/Documents/ETH/Master/Master_Thesis/LacI_TetR_Tup1/PMA1/',datestr(now, 'dd-mmm-yyyy'),'_rand_parameter_LacI_TetR_PMA1.mat'],'rand_parameter')

