
%Create the model
TetRmodel = sbiomodel('TetRmodel');

%Add a compartement
c = addcompartment(TetRmodel, 'comp');

%Add the system's species
s1 = addspecies(TetRmodel, 'TetR', 'InitialAmount', 0);
s2 = addspecies(TetRmodel, 'Citimmature', 'InitialAmount',0);
s3 = addspecies(TetRmodel, 'Citrine', 'InitialAmount',0);
s4 = addspecies(TetRmodel, 'TetRfree');

%Add parameters
p1 = addparameter(TetRmodel, 'k7tetTetR', 'Value', 50);
p2 = addparameter(TetRmodel, 'k7tetCit', 'Value', 50);
p3 = addparameter(TetRmodel, 'dTetR', 'Value', 0.5);
p4 = addparameter(TetRmodel, 'dCit', 'Value', 0.5);
p5 = addparameter(TetRmodel, 'kL7tet', 'Value', 1e-02);
p6 = addparameter(TetRmodel, 'thetaTetR', 'Value', 1);
p7 = addparameter(TetRmodel, 'KdTetR', 'Value', 0.5);
p8 = addparameter(TetRmodel, 'nTetR', 'Value', 2.5);
p9 = addparameter(TetRmodel, 'mu', 'Value', 0.0077);
p10 = addparameter(TetRmodel, 'kmaturation', 'Value', 0.0173);
p11 = addparameter(TetRmodel, 'atc', 'Value', 250);

%Add rate rules
raterule1 = addrule(TetRmodel, 'TetR = k7tetTetR*(kL7tet+(1-kL7tet)/(1+(TetRfree/thetaTetR)^nTetR)) - (dTetR+mu)*TetR', 'RuleType', 'rate');
raterule2 = addrule(TetRmodel, 'Citimmature = k7tetCit*(kL7tet+(1-kL7tet)/(1+(TetRfree/thetaTetR)^nTetR)) - (dCit+mu)*Citimmature', 'RuleType', 'rate');
raterule3 = addrule(TetRmodel, 'Citrine = kmaturation* Citimmature - (dCit +mu)*Citrine', 'RuleType', 'rate');

%Add repeated assignment rule
repassrule = addrule(TetRmodel, 'TetRfree = TetR/2 - KdTetR/2 - atc/2 + sqrt(KdTetR^2 + 2*KdTetR*TetR + 2*KdTetR*atc + TetR^2 - 2*TetR*atc + atc^2)/2', 'RuleType', 'repeatedAssignment');

%Simulate the model
[t,sd,species] = sbiosimulate(TetRmodel);

%Plot results
plot(t,sd);
legend(species);
xlabel('Time');
ylabel('Species Amount');