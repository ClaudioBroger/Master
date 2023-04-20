%Create model
LacI_TetR_Tup1_model = sbiomodel('LacI_TetR_Tup1_model');

%Add a compartment
c = addcompartment(LacI_TetR_Tup1_model, 'comp');

%Add the system's species
s1 = addspecies(LacI_TetR_Tup1_model, 'LacI', 'InitialAmount', 0);
s2 = addspecies(LacI_TetR_Tup1_model, 'Citrimmature', 'InitialAmount', 0);
s3 = addspecies(LacI_TetR_Tup1_model, 'Citrine', 'InitialAmount', 0);
s4 = addspecies(LacI_TetR_Tup1_model, 'TetRTup1', 'InitialAmount', 0);
s5 = addspecies(LacI_TetR_Tup1_model, 'LacIfree', 'Initialamount', 0);
s6 = addspecies(LacI_TetR_Tup1_model, 'Tup1free', 'InitialAmount', 0);
s7 = addspecies(LacI_TetR_Tup1_model, 'IPTG', 'InitialAmount', 0);
s8 = addspecies(LacI_TetR_Tup1_model, 'IPTG_InUnit', 'InitialAmount', 0);
s9 = addspecies(LacI_TetR_Tup1_model, 'aTc', 'InitialAmount', 0);
s10 = addspecies(LacI_TetR_Tup1_model, 'aTc_InUnit', 'InitialAmount',  0);
s11 = addspecies(LacI_TetR_Tup1_model, 'KdLacI_InUnit', 'InitialAmount', 0);
s12 = addspecies(LacI_TetR_Tup1_model, 'KdTetRTup1_InUnit', 'InitialAmount', 0);
s13 = addspecies(LacI_TetR_Tup1_model, 'TetR', 'InitialAmount', 0);
c14 = addspecies(LacI_TetR_Tup1_model, 'TetRfree', 'InitialAmount', 0);


%Add parameters
p1 = addparameter(LacI_TetR_Tup1_model, 'k7tetTetR', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p2 = addparameter(LacI_TetR_Tup1_model, 'kactTetR', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p3 = addparameter(LacI_TetR_Tup1_model, 'krnrTup1', 'Value', 0.0014, 'ValueUnits', '(molarity)/minute');
p4 = addparameter(LacI_TetR_Tup1_model, 'k7tetPMA1', 'Value', 0.0077, 'ValueUnits', '(molarity)/minute');
p5 = addparameter(LacI_TetR_Tup1_model, 'k7tetCit', 'Value', 0.3607, 'ValueUnits', '(molarity)/minute');
p6 = addparameter(LacI_TetR_Tup1_model, 'mu', 'Value', 0.0077, 'Valueunits', '1/minute');
p7 = addparameter(LacI_TetR_Tup1_model, 'kmaturation', 'Value', 0.0173, 'ValueUnits', '1/minute');
p8 = addparameter(LacI_TetR_Tup1_model, 'nMperUnit', 'Value', 1.6552, 'ValueUnits', 'dimensionless');
p9 = addparameter(LacI_TetR_Tup1_model, 'kactLacI', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p10 = addparameter(LacI_TetR_Tup1_model, 'k4lacLacI', 'Value', 0.0077, 'ValueUnits', '(molarity)/minute');
p11 = addparameter(LacI_TetR_Tup1_model, 'kactLacICit', 'Value', 0.4532, 'ValueUnits', '(molarity)/minute');
p12 = addparameter(LacI_TetR_Tup1_model, 'k3lacCit', 'Value', 0.3607, 'ValueUnits', '(molarity)/minute');
p13 = addparameter(LacI_TetR_Tup1_model, 'k4lacCit', 'Value', 0.3607, 'ValueUnits', '(molarity)/minute');
p14 = addparameter(LacI_TetR_Tup1_model, 'k3lacTup1', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p15 = addparameter(LacI_TetR_Tup1_model, 'dTetR', 'Value', 0.0014, 'ValueUnits', '1/minute');
p16 = addparameter(LacI_TetR_Tup1_model, 'dTup1', 'Value', 1.001, 'ValueUnits', '1/minute');
p17 = addparameter(LacI_TetR_Tup1_model, 'dTup1adh1', 'Value', 0.4532, 'ValueUnits', '1/minute');
p18 = addparameter(LacI_TetR_Tup1_model, 'dPMA1', 'Value', 0.3607, 'ValueUnits', '1/minute');
p19 = addparameter(LacI_TetR_Tup1_model, 'dLacI', 'Value', 0.3607, 'ValueUnits', '1/minute');
p20 = addparameter(LacI_TetR_Tup1_model, 'dLacICit', 'Value', 0.3607, 'ValueUnits', '1/minute');
p21 = addparameter(LacI_TetR_Tup1_model, 'kL7tet', 'Value', 0.0014, 'ValueUnits', 'dimensionless');
p22 = addparameter(LacI_TetR_Tup1_model, 'kl4lac', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p23 = addparameter(LacI_TetR_Tup1_model, 'kl3lac', 'Value', 1.0803, 'ValueUnits', 'dimensionless');
p24 = addparameter(LacI_TetR_Tup1_model, 'thetaTetR', 'Value', 0.0014, 'ValueUnits', 'molarity');
p25 = addparameter(LacI_TetR_Tup1_model, 'thetaTup1', 'Value', 1, 'ValueUnits', 'molarity');
p26 = addparameter(LacI_TetR_Tup1_model, 'thetaLacI', 'Value', 1, 'ValueUnits', 'molarity');
p27 = addparameter(LacI_TetR_Tup1_model, 'thetaLacI_2', 'Value', 1, 'ValueUnits', 'molarity');
p28 = addparameter(LacI_TetR_Tup1_model, 'thetaLacImut', 'Value', 1, 'ValueUnits', 'molarity');
p29 = addparameter(LacI_TetR_Tup1_model, 'thetaLacI3mut', 'Value', 1, 'ValueUnits', 'molarity');
p30 = addparameter(LacI_TetR_Tup1_model, 'thetaLacI3mut_2', 'Value', 1, 'ValueUnits', 'molarity');
p31 = addparameter(LacI_TetR_Tup1_model, 'KdTetR', 'Value', 0.5, 'ValueUnits', 'mole/liter');
p32 = addparameter(LacI_TetR_Tup1_model, 'KdLacI', 'Value', 0.5, 'ValueUnits', 'mole/liter');
p33 = addparameter(LacI_TetR_Tup1_model, 'nTetR', 'Value', 2.5, 'ValueUnits', 'dimensionless');
p34 = addparameter(LacI_TetR_Tup1_model, 'nTup1', 'Value', 2.5, 'ValueUnits', 'dimensionless');
p35 = addparameter(LacI_TetR_Tup1_model, 'nLacI', 'Value', 2, 'ValueUnits', 'dimensionless');
p36 = addparameter(LacI_TetR_Tup1_model, 'nLacI2', 'Value', 2, 'ValueUnits', 'dimensionless');
p37 = addparameter(LacI_TetR_Tup1_model, 'a', 'Value', 1, 'ValueUnits', 'dimensionless');
p38 = addparameter(LacI_TetR_Tup1_model, 'b', 'Value', 1, 'ValueUnits', 'dimensionless');
p39 = addparameter(LacI_TetR_Tup1_model, 'f', 'Value', 1, 'ValueUnits', 'dimensionless');
p40 = addparameter(LacI_TetR_Tup1_model, 'g', 'Value', 1, 'ValueUnits', 'dimensionless');
p41 = addparameter(LacI_TetR_Tup1_model, 'p1', 'Value', 5, 'ValueUnits', 'dimensionless');
p42 = addparameter(LacI_TetR_Tup1_model, 'p2', 'Value', 5, 'ValueUnits', 'dimensionless');


%scaling factor
scaling1 = addrule(LacI_TetR_Tup1_model, 'KdLacI_InUnit = KdLacI/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling2 = addrule(LacI_TetR_Tup1_model, 'IPTG_InUnit = IPTG/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling3 = addrule(LacI_TetR_Tup1_model, 'KdTetR_Tup1_InUnit = KdTetR_Tup1/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling4 = addrule(LacI_TetR_Tup1_model, 'ATC_InUnit = ATC/nMperUnit', 'RuleType', 'repeatedAssignment');

%raterules
raterule1 = addrule(LacI_TetR_Tup1_model, 'TetR = kactTetR - (dTetR + mu) * TetR', 'RuleType', 'rate');
raterule2 = addrule(LacI_TetR_Tup1_model, 'TetRTup1 = krnrTup1 + k3lacTup1 * (kl3lac + (1 - kl3lac)/(1 + (LacIfree/(thetaLacI_2 + thetalacI3mut_2 * b))^nLacI2)) - (dTup1 + dTup1adh1 + mu) * TetRTup1', 'RuleType', 'rate');
raterule3 = addrule(LacI_TetR_Tup1_model, 'Citrimmature = k7tetCit * (kL7tet + (1 - kL7tet)/(1 + (TetRfree * f / thetaTetR)^nTetR + (Tup1free * g / thetaTup1)^nTup1)) - (dCit + mu + kmaturation) * Citrimmature', 'RuleType', 'rate');
raterule4 = addrule(LacI_TetR_Tup1_model, 'Citrine = (kmaturation * Citrimmature) - (dCit + mu) * Citrine', 'RuleType', 'rate');
raterule5 = addrule(LacI_TetR_Tup1_model, 'LacI = (kactLacI + kactLacICit) - (dLacI + dLacICit + mu) * LacI', 'RuleType', 'rate');

%repeated assignment rules
repassrule1 = addrule(LacI_TetR_Tup1_model, 'LacIfree = max(0, LacI/2 - KdLacI_InUnit/2 - IPTG_InUnit/2 + sqrt(KdLacI_InUnit^2 + 2*KdLacI_InUnit*LacI + 2*KdLacI_InUnit*IPTG_InUnit + LacI^2 - 2*LacI*IPTG_InUnit + IPTG_InUnit^2)/2)', 'RuleType', 'repeatedAssignment');
repassrule2 = addrule(LacI_TetR_Tup1_model, 'TetRfree = max(TetR - (TetR*(TetR*TetRTup1 + TetRTup1*aTc_InUnit - TetRTup1*(KdTetRTup1_InUnit^2 + 2*KdTetRTup1_InUnit*TetR + 2*KdTetRTup1_InUnit*TetRTup1 + 2*KdTetRTup1_InUnit*aTc_InUnit + TetR^2 + 2*TetR*TetRTup1 - 2*TetR*aTc_InUnit + TetRTup1^2 - 2*TetRTup1*aTc_InUnit + aTc_InUnit^2)^(1/2) + TetRTup1^2 + KdTetRTup1_InUnit*TetRTup1))/(2*TetRTup1*(TetR + TetRTup1)),0);', 'RuleType', 'repeatedAssignment');
repassrule3 = addrule(LacI_TetR_Tup1_model, 'Tup1free = max(TetRTup1 - (TetR*TetRTup1 + TetRTup1*aTc_InUnit - TetRTup1*(KdTetRTup1_InUnit^2 + 2*KdTetRTup1_InUnit*TetR + 2*KdTetRTup1_InUnit*TetRTup1 + 2*KdTetRTup1_InUnit*aTc_InUnit + TetR^2 + 2*TetR*TetRTup1 - 2*TetR*aTc_InUnit + TetRTup1^2 - 2*TetRTup1*ATC_InUnit + aTc_InUnit^2)^(1/2) + TetRTup1^2 + KdTetRTup1_InUnit*TetRTup1)/(2*(TetR + TetRTup1)),0);', 'RuleType', 'repeatedAssignment');




