%% Create the model
%Create model object 
TetR = sbiomodel('TetR');

%Adding a compartement to the model object with the units in liters

compObj = addcompartment(TetR,'comp');
compObj.CapacityUnits = 'liter';

%Adding the models species

s = addspecies(TetR, 'TetRfree', 'InitialAmount', 0);
s1 = addspecies(TetR, 'TetR', 'InitialAmount', 0);
s2 = addspecies(TetR, 'Citimmature', 'InitialAmount', 0);
s3 = addspecies(TetR, 'Citrine', 'InitialAmount', 0);

%Adding model parameters

k7tetTetR = addparameter(TetR, 'k7tetTetR', 'Value', 50);
k7tetCit = addparameter(TetR, 'k7tetCit', 'Value', 50);
dTetR = addparameter(TetR, 'dTetR', 'Value', 0.5);
dCit = addparameter(TetR, 'dCit', 'Value', 0.5);
kL7tet = addparameter(TetR, 'kL7tet', 'Value', 1e-02);
thetaTetR = addparameter(TetR, 'thetaTetR', 'Value', 1);
KdTetR = addparameter(TetR, 'KdTetR', 'Value', 0.5);
nTetR = addparameter(TetR, 'nTetR', 'Value', 2.5);

%Adding Induction parameters

actAdded = addparameter(TetR, 'atcAdded', 'Value', 0);
indTime = addparameter(TetR, 'indTime', 'Value', 2000);

%Adding switching parameters

mu = addparameter(TetR, 'mu', 'Value', 0.0077);

%Mn to unit parameters

nMperUnit = addparameter(TetR, 'nMperUnit', 'Value', 1);

%Maturation time Citrine

kmaturation = addparameter(TetR, 'kmaturation', 'Value', 0.0173);

%Inducer concentration

atc = addparameter(TetR, 'atc', 'Value', 250);

%Converting to unit

KdTetR_InUnit = addparameter(TetR, 'KdTetR_InUnit', KdTetR/nMperUnit);

%Equation (rule) for TetR

r = addrule(TetR, 'TetRfree = ((TetR/2)-(KdTetR_InUnit/2)-(atc/2)+(KdTetR_InUnit^2+2*KdTetR_InUnit+2*KdTetR_InUnit*atc+TetR^2-2*TetR*atc+atc^2)^(1/2))/2', 'Ruletype', 'rate')

[t, sd, species] = sbiosimulate(TetR);
plot(t,sd);
legend(species);
xlabel('Time');
ylabel('Species Amount');
