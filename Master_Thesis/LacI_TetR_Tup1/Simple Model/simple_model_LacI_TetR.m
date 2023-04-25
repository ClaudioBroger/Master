% Simple model with three species

simple_model = sbiomodel('simple_model');

c = addcompartment(simple_model, 'comp');

s1 = addspecies(simple_model, 'LacI', 'InitialAmount', 0);
s2 = addspecies(simple_model, 'LacIfree', 'InitialAmount', 0);
s3 = addspecies(simple_model, 'TetRTup1', 'InitialAmount', 0);
s4 = addspecies(simple_model, 'TetRTup1free', 'InitialAmount', 0);
s5 = addspecies(simple_model, 'Citrimmature', 'InitialAmount', 0);
s6 = addspecies(simple_model, 'Citrine', 'InitialAmount', 0);
s7 = addspecies(simple_model, 'IPTG', 'InitialAmount', 0);
s8 = addspecies(simple_model, 'IPTG_InUnit', 'InitialAmount', 0);
s9 = addspecies(simple_model, 'KdLacI_InUnit', 'InitialAmount', 0);
s10 = addspecies(simple_model, 'KdTetR_InUnit', 'InitialAmount', 0);

p1 = addparameter(simple_model, 'kLacI', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p2 = addparameter(simple_model, 'dLacI', 'Value', 0.0014, 'ValueUnits', '1/minute');
p3 = addparameter(simple_model, 'kTetRTup1', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p4 = addparameter(simple_model, 'dTetRTup1', 'Value', 0.0014, 'ValueUnits', '1/minute');
p5 = addparameter(simple_model, 'degtag', 'Value', 0.0077, 'ValueUnits', '1/minute');
p6 = addparameter(simple_model, 'LacIrep', 'Value', 0.3607, 'ValueUnits', 'molarity');
p7 = addparameter(simple_model, 'LacIrep2', 'Value', 0.3607, 'ValueUnits', 'molarity');
p8 = addparameter(simple_model, 'LacIrep3', 'Value', 0.3607, 'ValueUnits', 'molarity');
p9 = addparameter(simple_model, 'kCit', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p10 = addparameter(simple_model, 'dCit', 'Value', 0.0077, 'ValueUnits', '1/minute');
p11 = addparameter(simple_model, 'CitL', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p12 = addparameter(simple_model, 'TetRTup1L', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p13 = addparameter(simple_model, 'nLacI', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p14 = addparameter(simple_model, 'nTetRTup1', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p15 = addparameter(simple_model, 'KdLacI', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p16 = addparameter(simple_model, 'KdTetR', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p17 = addparameter(simple_model, 'kmaturation', 'Value', 0.0173, 'ValueUnits', '1/minute');
p18 = addparameter(simple_model, 'nMperUnit', 'Value', 1.6552, 'ValueUnits', 'dimensionless');
p19 = addparameter(simple_model, 'mu', 'Value', 0.0077, 'Valueunits', '1/minute');
p20 = addparameter(simple_model, 'TetRTup1rep', 'Value', 0.3607, 'ValueUnits', 'molarity');

scaling1 = addrule(simple_model, 'KdLacI_InUnit = KdLacI/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling2 = addrule(simple_model, 'IPTG_InUnit = IPTG/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling3 = addrule(simple_model, 'KdTetR_InUnit = KdTetR/nMperUnit', 'RuleType', 'repeatedAssignment');

raterule1 = addrule(simple_model, 'LacI = kLacI - (dLacI + mu) * LacI', 'RuleType', 'rate');
raterule2 = addrule(simple_model, 'TetRTup1 = kTetRTup1 * (TetRTup1L + (1 - TetRTup1L)/(1 + (LacIfree/(LacIrep+LacIrep2+LacIrep3))^nLacI)) - (dTetRTup1 + mu + degtag) * TetRTup1', 'RuleType', 'rate');
raterule3 = addrule(simple_model, 'Citrimmature = kCit * (CitL + (1 - CitL)/(1 + (TetRTup1/TetRTup1rep)^nTetRTup1)) - (dCit + mu + kmaturation) * Citrimmature', 'RuleType', 'rate');
raterule4 = addrule(simple_model, 'Citrine = kmaturation * Citrimmature - (dCit + mu) * Citrine', 'RuleType', 'rate');

repassrule1 = addrule(simple_model, 'LacIfree = max(0, LacI/2 - KdLacI_InUnit/2 - IPTG_InUnit/2 + sqrt(KdLacI_InUnit^2 + 2*KdLacI_InUnit*LacI + 2*KdLacI_InUnit*IPTG_InUnit + LacI^2 - 2*LacI*IPTG_InUnit + IPTG_InUnit^2)/2)', 'RuleType', 'repeatedAssignment');




