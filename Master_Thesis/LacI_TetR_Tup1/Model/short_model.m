%%Short LacI TetRTup1 Cit model

short_LacI_TetR_model = sbiomodel('short_LacI_TetR_model');

c = addcompartment(short_LacI_TetR_model, 'comp');

s1 = addspecies(short_LacI_TetR_model, 'LacI', 'InitialAmount', 0);
s2 = addspecies(short_LacI_TetR_model, 'TetRTup1', 'InitialAmount', 0);
s3 = addspecies(short_LacI_TetR_model, 'LacIfree', 'InitialAmount', 0);
s4 = addspecies(short_LacI_TetR_model, 'TetRTup1free', 'InitialAmount', 0);
s5 = addspecies(short_LacI_TetR_model, 'Citrimmature', 'InitialAmount', 0);
s6 = addspecies(short_LacI_TetR_model, 'Citrine', 'InitialAmount', 0);
s7 = addspecies(short_LacI_TetR_model, 'IPTG', 'InitialAmount', 0);
s8 = addspecies(short_LacI_TetR_model, 'IPTG_InUnit', 'InitialAmount', 0);
s9 = addspecies(short_LacI_TetR_model, 'aTc', 'InitialAmount', 0);
s10 = addspecies(short_LacI_TetR_model, 'aTc_InUnit', 'InitialAmount', 0);
s11 = addspecies(short_LacI_TetR_model, 'KdLacI_InUnit', 'InitialAmount', 0);
s12 = addspecies(short_LacI_TetR_model, 'KdTetRTup1_InUnit', 'InitialAmount', 0);

p1 = addparameter(short_LacI_TetR_model, 'kLacI', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p2 = addparameter(short_LacI_TetR_model, 'kTetRTup1', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p3 = addparameter(short_LacI_TetR_model, 'kCit', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p4 = addparameter(short_LacI_TetR_model, 'dLacI', 'Value', 0.0014, 'ValueUnits', '1/minute');
p5 = addparameter(short_LacI_TetR_model, 'dCit', 'Value', 0.0077, 'ValueUnits', '1/minute');
p6 = addparameter(short_LacI_TetR_model, 'dTetRTup1', 0.0014, 'ValueUnits', '1/minute');
p7 = addparameter(short_LacI_TetR_model, 'LacIrep', 'Value', 0.3607, 'ValueUnits', 'molarity');
p8 = addparameter(short_LacI_TetR_model, 'mu', 'Value', 0.0077, 'Valueunits', '1/minute');
p9 = addparameter(short_LacI_TetR_model, 'kmaturation', 'Value', 0.0173, 'ValueUnits', '1/minute');
p10 = addparameter(short_LacI_TetR_model, 'nMperUnit', 'Value', 1.6552, 'ValueUnits', 'dimensionless');
p11 = addparameter(short_LacI_TetR_model, 'nLacI', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p12 = addparameter(short_LacI_TetR_model, 'CitL', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p13 = addparameter(short_LacI_TetR_model, 'KdLacI', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p14 = addparameter(short_LacI_TetR_model, 'LacIrep2', 'Value', 0.3607, 'ValueUnits', 'molarity');
p15 = addparameter(short_LacI_TetR_model, 'LacIrep3', 'Value', 0.3607, 'ValueUnits', 'molarity');
p16 = addparameter(short_LacI_TetR_model, 'nTetRTup1', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p17 = addparameter(short_LacI_TetR_model, 'TetRTup1_L', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p18 = addparameter(short_LacI_TetR_model, 'KdTetRTup1', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p19 = addparameter(short_LacI_TetR_model, 'TetRTup1rep', 'Value', 0.3607, 'ValueUnits', 'molarity');
p20 = addparameter(short_LacI_TetR_model, 'degtag', 'Value', 0.0014, 'ValueUnits', '1/minute');

scaling1 = addrule(short_LacI_TetR_model, 'KdLacI_InUnit = KdLacI/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling2 = addrule(short_LacI_TetR_model, 'IPTG_InUnit = IPTG/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling3 = addrule(short_LacI_TetR_model, 'KdTetRTup1_InUnit = KdTetRTup1/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling4 = addrule(short_LacI_TetR_model, 'aTc_InUnit = aTc/nMperUnit', 'RuleType', 'repeatedAssignment');

raterule1 = addrule(short_LacI_TetR_model, 'LacI = kLacI - (dLacI + mu) * LacI', 'RuleType', 'rate');
raterule2 = addrule(short_LacI_TetR_model, 'Citrimmature = kCit * (CitL + (1 - CitL)/(1 + (TetRTup1free/(TetRTup1rep))^nTetRTup1)) - (dCit + mu + kmaturation) * Citrimmature', 'RuleType', 'rate');
raterule3 = addrule(short_LacI_TetR_model, 'Citrine = kmaturation * Citrimmature - (dCit + mu) * Citrine', 'RuleType', 'rate');
raterule4 = addrule(short_LacI_TetR_model, 'TetRTup1 = kTetRTup1 * (TetRTup1_L + (1 - TetRTup1_L)/(1 + (LacIfree/(LacIrep + LacIrep2 + LacIrep3))^nLacI)) - (dTetRTup1 + mu + degtag) * TetRTup1', 'RuleType', 'rate');

repassrule1 = addrule(short_LacI_TetR_model, 'LacIfree = max(0, LacI/2 - KdLacI_InUnit/2 - IPTG_InUnit/2 + sqrt(KdLacI_InUnit^2 + 2*KdLacI_InUnit*LacI + 2*KdLacI_InUnit*IPTG_InUnit + LacI^2 - 2*LacI*IPTG_InUnit + IPTG_InUnit^2)/2)', 'RuleType', 'repeatedAssignment');
repassrule2 = addrule(short_LacI_TetR_model, 'TetRTup1 = max(0, TetRTup1/2 - KdTetRTup1_InUnit/2 - aTc_InUnit/2 + sqrt(KdTetRTup1_InUnit^2 + 2*KdTetRTup1_InUnit*TetRTup1 + 2*KdTetRTup1_InUnit*aTc_InUnit + TetRTup1^2 - 2*TetRTup1*aTc_InUnit + aTc_InUnit^2)/2)', 'RuleType', 'repeatedAssignment');





