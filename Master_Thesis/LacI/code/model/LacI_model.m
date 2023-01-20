%Create the model
LacImodel = sbiomodel('LacImodel');

%Add a compartement
c = addcompartment(LacImodel, 'comp');

%Add the system's species
s1 = addspecies(LacImodel, 'LacI', 'InitialAmount', 0);
s2 = addspecies(LacImodel, 'Citimmature', 'InitialAmount', 0);
s3 = addspecies(LacImodel, 'Citimmature_1', 'InitialAmount', 0);
s4 = addspecies(LacImodel, 'Citrine', 'InitialAmount', 0);
s5 = addspecies(LacImodel, 'LacIfree');
s6 = addspecies(LacImodel, 'IPTG', 'InitialAmount', 0);
s7 = addspecies(LacImodel, 'IPTG_InUnit', 'InitialAmount', 0);
s8 = addspecies(LacImodel, 'KdLacI_InUnit', 'InitialAmount', 0);

%Add parameters
p1 = addparameter(LacImodel, 'PAct1_LacI', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p2 = addparameter(LacImodel, 'P4Lacn_cit', 'Value', 1.3389, 'ValueUnits', '(molarity)/minute');
p3 = addparameter(LacImodel, 'dLacI', 'Value', 0.0014, 'ValueUnits', '1/minute');
p4 = addparameter(LacImodel, 'dCit', 'Value', 0.0077, 'ValueUnits', '1/minute');
p5 = addparameter(LacImodel, 'PAct1_LacI_L','Value', 0.0077, 'ValueUnits', 'dimensionless');
p6 = addparameter(LacImodel, 'LacI_rep_Cit', 'Value', 0.3607, 'ValueUnits', 'molarity');
p7 = addparameter(LacImodel, 'KdLacI', 'Value', 0.4532, 'ValueUnits', 'mole/liter');
p8 = addparameter(LacImodel, 'nLacI', 'Value', 1.001, 'ValueUnits', 'dimensionless');
p9 = addparameter(LacImodel, 'mu', 'Value', 0.0077, 'Valueunits', '1/minute');
p10 = addparameter(LacImodel, 'kmaturation', 'Value', 0.0173, 'ValueUnits', '1/minute');
p11 = addparameter(LacImodel, 'nMperUnit', 'Value', 1.6552, 'ValueUnits', 'dimensionless');
p12 = addparameter(LacImodel, 'LacI_rep_Cit_W220F', 'Value', 3.607, 'ValueUnits', 'molarity');
p13 = addparameter(LacImodel, 'P_4Lacn_LacI', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p14 = addparameter(LacImodel, 'P_4Lacn_LacI_L', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p15 = addparameter(LacImodel, 'LacI_rep', 'Value', 0.3607, 'ValueUnits', 'molarity');
p16 = addparameter(LacImodel, 'Silence_LacI_rep', 'Value', 1, 'ValueUnits', 'dimensionless');
p17 = addparameter(LacImodel, 'pt7_LacI', 'Value', 1.0803, 'ValueUnits', '(molarity)/minute');
p18 = addparameter(LacImodel, 'P3_Lacn_5_cit', 1.3389, 'ValueUnits', '(molarity)/minute');
p19 = addparameter(LacImodel, 'P3_Lacn_5_cit_L', 'Value', 0.0077, 'ValueUnits', 'dimensionless');
p20 = addparameter(LacImodel, 'dLacI_pt7', 'Value', 0.0014, 'ValueUnits', '1/minute');
p21 = addparameter(LacImodel, 'nLacI_P3', 'Value', 1.001, 'ValueUnits', 'dimensionless');

%Add scaling factor to KdLacI and IPTG
scaling1 = addrule(LacImodel, 'KdLacI_InUnit = KdLacI/nMperUnit', 'RuleType', 'repeatedAssignment');
scaling2 = addrule(LacImodel, 'IPTG_InUnit = IPTG/nMperUnit', 'RuleType', 'repeatedAssignment');

%Add rate rules
raterule1 = addrule(LacImodel, 'LacI = PAct1_LacI + pt7_LacI + P_4Lacn_LacI*(P_4Lacn_LacI_L+PAct1_LacI_L + (1- P_4Lacn_LacI_L+PAct1_LacI_L)/(1 + (LacIfree/LacI_rep)^(nLacI+nLacI_P3)))-(dLacI+dLacI_pt7+mu)*LacI', 'RuleType', 'rate');
raterule2 = addrule(LacImodel, 'Citimmature = P4Lacn_cit*(P_4Lacn_LacI_L+(1- P_4Lacn_LacI_L)/(1+(LacIfree/(LacI_rep_Cit+LacI_rep_Cit_W220F+(Silence_LacI_rep*LacI_rep)))^nLacI)) - (dCit+mu)*Citimmature', 'RuleType', 'rate');
raterule3 = addrule(LacImodel, 'Citimmature_1 = P3_Lacn_5_cit*(P3_Lacn_5_cit_L+(1-P3_Lacn_5_cit_L)/(1+(LacIfree/(LacI_rep_Cit+LacI_rep_Cit_W220F+(Silence_LacI_rep*LacI_rep)))^nLacI_P3)) - (dCit+mu)*Citimmature', 'RuleType', 'rate');
raterule4 = addrule(LacImodel, 'Citrine = kmaturation*(Citimmature+Citimmature_1) - (dCit+mu)*Citrine', 'RuleType', 'rate');

%Add repeated assignment rule
repassrule = addrule(LacImodel, 'LacIfree = max(0, LacI/2 - KdLacI_InUnit/2 - IPTG_InUnit/2 + sqrt(KdLacI_InUnit^2 + 2*KdLacI_InUnit*LacI + 2*KdLacI_InUnit*IPTG_InUnit + LacI^2 - 2*LacI*IPTG_InUnit + IPTG_InUnit^2)/2)', 'RuleType', 'repeatedAssignment');