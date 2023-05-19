% Simple model with three species

LacI_TetRTup1_model = sbiomodel('LacI_TetRTup1_model');

c = addcompartment(LacI_TetRTup1_model, 'comp');

s1 = addspecies(LacI_TetRTup1_model, 'LacI', 'InitialAmount', 0);
s2 = addspecies(LacI_TetRTup1_model, 'LacIfree', 'InitialAmount', 0);
s3 = addspecies(LacI_TetRTup1_model, 'TetRTup1', 'InitialAmount', 0);
s4 = addspecies(LacI_TetRTup1_model, 'TetRTup1free', 'InitialAmount', 0);
s5 = addspecies(LacI_TetRTup1_model, 'Citrimmature', 'InitialAmount', 0);
s6 = addspecies(LacI_TetRTup1_model, 'Citrine', 'InitialAmount', 0);
s7 = addspecies(LacI_TetRTup1_model, 'IPTG', 'InitialAmount', 0);
s8 = addspecies(LacI_TetRTup1_model, 'IPTG_InUnit', 'InitialAmount', 0);
s9 = addspecies(LacI_TetRTup1_model, 'KdLacI_InUnit', 'InitialAmount', 0);
s10 = addspecies(LacI_TetRTup1_model, 'KdTetR_InUnit', 'InitialAmount', 0);
s11 = addspecies(LacI_TetRTup1_model, 'aTc', 'InitialAmount', 0);
s12 = addspecies(LacI_TetRTup1_model, 'aTc_InUnit', 'InitialAmount', 0);
s13 = addspecies(LacI_TetRTup1_model, 'TetR', 'InitialAmount', 0);

p1 = addparameter(LacI_TetRTup1_model, 'kLacI', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p2 = addparameter(LacI_TetRTup1_model, 'dLacI', 'Value', 0.0014, 'ValueUnits', '1/minute');
p3 = addparameter(LacI_TetRTup1_model, 'kTetRTup1', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p4 = addparameter(LacI_TetRTup1_model, 'dTetRTup1', 'Value', 0.0014, 'ValueUnits', '1/minute');
p5 = addparameter(LacI_TetRTup1_model, 'degtag', 'Value', 0.0077, 'ValueUnits', '1/minute');
p6 = addparameter(LacI_TetRTup1_model, 'LacIrep', 'Value', 0.3607, 'ValueUnits', 'molarity');
p7 = addparameter(LacI_TetRTup1_model, 'LacIrep2', 'Value', 0.3607, 'ValueUnits', 'molarity');
p8 = addparameter(LacI_TetRTup1_model, 'LacIrep3', 'Value', 0.3607, 'ValueUnits', 'molarity');
p9 = addparameter(LacI_TetRTup1_model, 'kCit', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p10 = addparameter(LacI_TetRTup1_model, 'dCit', 'Value', 0.0077, 'ValueUnits', '1/minute');
p11 = addparameter(LacI_TetRTup1_model, 'CitL', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p12 = addparameter(LacI_TetRTup1_model, 'TetRTup1L', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p13 = addparameter(LacI_TetRTup1_model, 'nLacI', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p14 = addparameter(LacI_TetRTup1_model, 'nTetRTup1', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p15 = addparameter(LacI_TetRTup1_model, 'KdLacI', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p16 = addparameter(LacI_TetRTup1_model, 'KdTetR', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p17 = addparameter(LacI_TetRTup1_model, 'kmaturation', 'Value', 0.0173, 'ValueUnits', '1/minute');
p18 = addparameter(LacI_TetRTup1_model, 'nMperUnit', 'Value', 1.6552, 'ValueUnits', 'dimensionless');
p19 = addparameter(LacI_TetRTup1_model, 'mu', 'Value', 0.0077, 'Valueunits', '1/minute');
p20 = addparameter(LacI_TetRTup1_model, 'TetRTup1rep', 'Value', 0.3607, 'ValueUnits', 'molarity');
p21 = addparameter(LacI_TetRTup1_model, 'kTetR', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p22 = addparameter(LacI_TetRTup1_model, 'dTetR', 'Value', 0.0014, 'ValueUnits', '1/minute');

scaling1 = addrule(LacI_TetRTup1_model, 'KdLacI_InUnit = KdLacI/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling2 = addrule(LacI_TetRTup1_model, 'IPTG_InUnit = IPTG/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling3 = addrule(LacI_TetRTup1_model, 'KdTetR_InUnit = KdTetR/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling4 = addrule(LacI_TetRTup1_model, 'aTc_InUnit = aTc/nMperUnit', 'RuleType', 'repeatedAssignment');

raterule1 = addrule(LacI_TetRTup1_model, 'LacI = kLacI - (dLacI + mu) * LacI', 'RuleType', 'rate');
raterule2 = addrule(LacI_TetRTup1_model, 'TetRTup1 = kTetRTup1 * (TetRTup1L + (1 - TetRTup1L)/(1 + (LacIfree/(LacIrep+LacIrep2+LacIrep3))^nLacI)) - (dTetRTup1 + mu + degtag) * TetRTup1', 'RuleType', 'rate');
raterule3 = addrule(LacI_TetRTup1_model, 'Citrimmature = kCit * (CitL + (1 - CitL)/(1 + (TetRTup1free/TetRTup1rep)^nTetRTup1)) - (dCit + mu + kmaturation) * Citrimmature', 'RuleType', 'rate');
raterule4 = addrule(LacI_TetRTup1_model, 'Citrine = kmaturation * Citrimmature - (dCit + mu) * Citrine', 'RuleType', 'rate');
raterule5 = addrule(LacI_TetRTup1_model, 'TetR = kTetR - (dTetR + mu) * TetR', 'RuleType', 'rate');

repassrule1 = addrule(LacI_TetRTup1_model, 'LacIfree = max(0, LacI/2 - KdLacI_InUnit/2 - IPTG_InUnit/2 + (KdLacI_InUnit^2 + 2*KdLacI_InUnit*LacI + 2*KdLacI_InUnit*IPTG_InUnit + LacI^2 - 2*LacI*IPTG_InUnit + IPTG_InUnit^2)^(1/2)/2)', 'RuleType', 'repeatedAssignment');
repassrule2 = addrule(LacI_TetRTup1_model, 'TetRTup1free = max(TetRTup1 - (TetR*TetRTup1 + TetRTup1*aTc_InUnit - TetRTup1*(KdTetR_InUnit^2 + 2*KdTetR_InUnit*TetR + 2*KdTetR_InUnit*TetRTup1 + 2*KdTetR_InUnit*aTc_InUnit + TetR^2 + 2*TetR*TetRTup1 - 2*TetR*aTc_InUnit + TetRTup1^2 - 2*TetRTup1*aTc_InUnit + aTc_InUnit^2)^(1/2) + TetRTup1^2 + KdTetR_InUnit*TetRTup1)/(2*(TetR + TetRTup1)),0)', 'RuleType', 'repeatedAssignment');
%repassrule2 = addrule(LacI_TetRTup1_model, 'TetRTup1free = max(0, TetRTup1/2 - KdTetR_InUnit/2 - aTc_InUnit/2 + (KdTetR_InUnit^2 + 2*KdTetR_InUnit*TetRTup1 + 2*KdTetR_InUnit*aTc_InUnit + TetRTup1^2 - 2*TetRTup1*aTc_InUnit + aTc_InUnit^2)^(1/2)/2)', 'RuleType', 'repeatedAssignment');



