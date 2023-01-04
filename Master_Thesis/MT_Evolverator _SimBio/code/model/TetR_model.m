%% Example of how to add dimensions

%Create the model
TetRmodel = sbiomodel('TetRmodel');

%Add a compartement
c = addcompartment(TetRmodel, 'comp');

%Add the system's species
s1 = addspecies(TetRmodel, 'TetR', 'InitialAmount', 0);
s2 = addspecies(TetRmodel, 'Citimmature', 'InitialAmount',0);
s3 = addspecies(TetRmodel, 'Citrine', 'InitialAmount',0);
s4 = addspecies(TetRmodel, 'TetRfree');
s5 = addspecies(TetRmodel,'atc', 'InitialAmount', 0);
s6 = addspecies(TetRmodel, 'atc_InUnit', 'InitialAmount', 0);
s7 = addspecies(TetRmodel, 'KdTetR_InUnit', 'InitialAmount', 0);

%Add parameters
p1 = addparameter(TetRmodel, 'k7tetTetR', 'Value', 1.0803, 'ValueUnits','(molarity)/minute');
p2 = addparameter(TetRmodel, 'k7tetCit', 'Value', 1.3389, 'ValueUnits','(molarity)/minute');
p3 = addparameter(TetRmodel, 'dTetR', 'Value', 0.0014, 'ValueUnits', '1/minute');
p4 = addparameter(TetRmodel, 'dCit', 'Value', 0.0077, 'ValueUnits', '1/minute');
p5 = addparameter(TetRmodel, 'kL7tet', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p6 = addparameter(TetRmodel, 'thetaTetR', 'Value', 0.3607, 'ValueUnits','molarity');
p7 = addparameter(TetRmodel, 'KdTetR', 'Value', 0.4532, 'ValueUnits','mole/liter');
p8 = addparameter(TetRmodel, 'nTetR', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p9 = addparameter(TetRmodel, 'mu', 'Value', 0.0077, 'ValueUnits', '1/minute');
p10 = addparameter(TetRmodel, 'kmaturation', 'Value', 0.0173, 'ValueUnits', '1/minute');
p11 = addparameter(TetRmodel, 'nMperUnit', 'Value', 1.6552, 'ValueUnits', 'dimensionless');

%Add scaling factor to KdTetR and atc
scaling1 = addrule(TetRmodel, 'KdTetR_InUnit = KdTetR/nMperUnit','RuleType', 'repeatedAssignment');
scaling2 = addrule(TetRmodel, 'atc_InUnit = atc/nMperUnit', 'RuleType', 'repeatedAssignment');

%Add rate rules
raterule1 = addrule(TetRmodel, 'TetR = k7tetTetR*(kL7tet+(1-kL7tet)/(1+(TetRfree/thetaTetR)^nTetR)) - (dTetR+mu)*TetR', 'RuleType', 'rate');
raterule2 = addrule(TetRmodel, 'Citimmature = k7tetCit*(kL7tet+(1-kL7tet)/(1+(TetRfree/thetaTetR)^nTetR)) - (dCit+mu)*Citimmature', 'RuleType', 'rate');
raterule3 = addrule(TetRmodel, 'Citrine = kmaturation* Citimmature - (dCit +mu)*Citrine', 'RuleType', 'rate');

%Add repeated assignment rule
%repassrule = addrule(TetRmodel, piecewise('TetR < 1e-04', 'TetRfree = 0', 'TetR >= 1e-04', 'TetRfree = TetR/2 - KdTetR_InUnit/2 - atc_InUnit/2 + sqrt(KdTetR_InUnit^2 + 2*KdTetR_InUnit*TetR + 2*KdTetR_InUnit*atc_InUnit + TetR^2 - 2*TetR*atc_InUnit + atc_InUnit^2)/2'), 'RuleType', 'repeatedAssignment');
%repassrule = addrule(TetRmodel, 'TetRfree = piecewise(TetR >= 0.0001, TetR/2 - KdTetR_InUnit/2 - atc_InUnit/2 + sqrt(KdTetR_InUnit^2 + 2*KdTetR_InUnit*TetR + 2*KdTetR_InUnit*atc_InUnit + TetR^2 - 2*TetR*atc_InUnit + atc_InUnit^2)/2, TetR < 0.0001, TetRfree = 0)', 'RuleType', 'repeatedAssignment');
%repassrule = addrule(TetRmodel, 'TetRfree = piecewise(TetR >= 0.0001, TetR/2 - KdTetR_InUnit/2 - atc_InUnit/2 + sqrt(KdTetR_InUnit^2 + 2*KdTetR_InUnit*TetR + 2*KdTetR_InUnit*atc_InUnit + TetR^2 - 2*TetR*atc_InUnit + atc_InUnit^2)/2, TetR < 0.0001, TetRfree = 0)', 'RuleType', 'repeatedAssignment');
repassrule = addrule(TetRmodel, 'TetRfree = max(0, TetR/2 - KdTetR_InUnit/2 - atc_InUnit/2 + sqrt(KdTetR_InUnit^2 + 2*KdTetR_InUnit*TetR + 2*KdTetR_InUnit*atc_InUnit + TetR^2 - 2*TetR*atc_InUnit + atc_InUnit^2)/2)', 'RuleType', 'repeatedAssignment');