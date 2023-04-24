%TetR_Tup1 model

TetR_Tup1_model = sbiomodel('TetR_Tup1_model');

c = addcompartment(TetR_Tup1_model, 'comp');

s1 = addspecies(TetR_Tup1_model, 'TetRTup1', 'InitialAmount', 0);
s2 = addspecies(TetR_Tup1_model, 'Citrimmature', 'InitialAmount', 0);
s3 = addspecies(TetR_Tup1_model, 'Citrine', 'InitialAmount', 0);
s4 = addspecies(TetR_Tup1_model, 'TetRfree', 'InitialAmount', 0);
s5 = addspecies(TetR_Tup1_model, 'Tup1free', 'InitialAmount', 0);
s6 = addspecies(TetR_Tup1_model, 'aTc', 'InitialAmount', 0);
s7 = addspecies(TetR_Tup1_model, 'aTc_InUnit', 'InitialAmount', 0);
s8 = addspecies(TetR_Tup1_model, 'KdTetRTup1_InUnit', 'InitialAmount', 0);

p1 = addparameter(TetR_Tup1_model, 'krnrTup1', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p2 = addparameter(TetR_Tup1_model, 'k7tetCit', 'Value',1.3389, 'ValueUnits', '(molarity)/minute');
p3 = addparameter(TetR_Tup1_model, 'dTup1', 'Value', 0.0014, 'ValueUnits', '1/minute');
p4 = addparameter(TetR_Tup1_model, 'dCit', 'Value', 0.0077, 'ValueUnits', '1/minute');
p5 = addparameter(TetR_Tup1_model, 'thetaTup1', 'Value', 0.3607, 'ValueUnits', 'molarity');
p6 = addparameter(TetR_Tup1_model, 'mu', 'Value', 0.0077, 'Valueunits', '1/minute');
p7 = addparameter(TetR_Tup1_model, 'kmaturation', 'Value', 0.0173, 'ValueUnits', '1/minute');
p8 = addparameter(TetR_Tup1_model, 'nMperUnit', 'Value', 1.6552, 'ValueUnits', 'dimensionless');
p9 = addparameter(TetR_Tup1_model, 'nTup1', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p10 = addparameter(TetR_Tup1_model, 'kL7tet', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p11 = addparameter(TetR_Tup1_model, 'KdTetR', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p12 = addparameter(TetR_Tup1_model, 'dTup1adh1', 'Value', 0.0077, 'ValueUnits', 'molarity');
p13 = addparameter(TetR_Tup1_model, 'kactTetR', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p14 = addparameter(TetR_Tup1_model, 'dTetR', 'Value', 0.0077, 'ValueUnits', '1/minute');
p15 = addparameter(TetR_Tup1_model, 'thetaTetR', 'Value', 0.3607, 'ValueUnits', 'molarity');
p16 = addparameter(TetR_Tup1_model, 'nTetR', 'Value', 1.001, 'ValueUnits', 'dimensionless');

scaling1 = addrule(TetR_Tup1_model, 'KdTetRTup1_InUnit = KdTetR/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling2 = addrule(TetR_Tup1_model, 'aTc_InUnit = aTc/nMperUnit', 'RuleType', 'repeatedAssignment');

raterule1 = addrule(TetR_Tup1_model, 'TetRTup1 = krnrTup1 - (dTup1 + mu + dTup1adh1) * TetRTup1', 'RuleType', 'rate');
raterule2 = addrule(TetR_Tup1_model,  'Citrimmature = k7tetCit * (kL7tet + (1 - kL7tet)/(1 + (TetRfree/ thetaTetR)^nTetR + (Tup1free/ thetaTup1)^nTup1)) - (dCit + mu + kmaturation) * Citrimmature', 'RuleType', 'rate');
raterule3 = addrule(TetR_Tup1_model, 'Citrine = (kmaturation * Citrimmature) - (dCit + mu) * Citrine', 'RuleType', 'rate');
raterule4 = addrule(TetR_Tup1_model, 'TetR = kactTetR - (dTetR + mu) * TetR', 'RuleType', 'rate');


repassrule1 = addrule(TetR_Tup1_model, 'TetRfree = max(TetR - (TetR*(TetR*TetRTup1 + TetRTup1*aTc_InUnit - TetRTup1*(KdTetRTup1_InUnit^2 + 2*KdTetRTup1_InUnit*TetR + 2*KdTetRTup1_InUnit*TetRTup1 + 2*KdTetRTup1_InUnit*aTc_InUnit + TetR^2 + 2*TetR*TetRTup1 - 2*TetR*aTc_InUnit + TetRTup1^2 - 2*TetRTup1*aTc_InUnit + aTc_InUnit^2)^(1/2) + TetRTup1^2 + KdTetRTup1_InUnit*TetRTup1))/(2*TetRTup1*(TetR + TetRTup1)),0);', 'RuleType', 'repeatedAssignment');
repassrule2 = addrule(TetR_Tup1_model, 'Tup1free = max(TetRTup1 - (TetR*TetRTup1 + TetRTup1*aTc_InUnit - TetRTup1*(KdTetRTup1_InUnit^2 + 2*KdTetRTup1_InUnit*TetR + 2*KdTetRTup1_InUnit*TetRTup1 + 2*KdTetRTup1_InUnit*aTc_InUnit + TetR^2 + 2*TetR*TetRTup1 - 2*TetR*aTc_InUnit + TetRTup1^2 - 2*TetRTup1*ATC_InUnit + aTc_InUnit^2)^(1/2) + TetRTup1^2 + KdTetRTup1_InUnit*TetRTup1)/(2*(TetR + TetRTup1)),0);', 'RuleType', 'repeatedAssignment');


