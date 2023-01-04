function [output] = odefile(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolverator with vectorization (test)
% Generated: 08-Mar-2021 15:03:11
% 
% [output] = odefile() => output = initial conditions in column vector
% [output] = odefile('states') => output = state names in cell-array
% [output] = odefile('algebraic') => output = algebraic variable names in cell-array
% [output] = odefile('parameters') => output = parameter names in cell-array
% [output] = odefile('parametervalues') => output = parameter values in column vector
% [output] = odefile('variablenames') => output = variable names in cell-array
% [output] = odefile('variableformulas') => output = variable formulas in cell-array
% [output] = odefile(time,statevector) => output = time derivatives in column vector
% 
% State names and ordering:
% 
% statevector(1): LexA1
% statevector(2): LexA2
% statevector(3): LexA3
% statevector(4): LexA4
% statevector(5): Mutator1
% statevector(6): Mutator2
% statevector(7): Mutator3
% statevector(8): Mutator4
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global time
parameterValuesNew = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
	% Return initial conditions of the state variables (and possibly algebraic variables)
	output = [0, 0, 0, 0, 0, 0, 0, 0];
	output = output(:);
	return
elseif nargin == 1,
	if strcmp(varargin{1},'states'),
		% Return state names in cell-array
		output = {'LexA1', 'LexA2', 'LexA3', 'LexA4', 'Mutator1', 'Mutator2', 'Mutator3', 'Mutator4'};
	elseif strcmp(varargin{1},'algebraic'),
		% Return algebraic variable names in cell-array
		output = {};
	elseif strcmp(varargin{1},'parameters'),
		% Return parameter names in cell-array
		output = {'kon1', 'kon2', 'kon3', 'kon4', 'N', 'p', 'd'};
	elseif strcmp(varargin{1},'parametervalues'),
		% Return parameter values in column vector
		output = [10, 100, 0.1, 1, 4, 1, 0.1];
	elseif strcmp(varargin{1},'variablenames'),
		% Return variable names in cell-array
		output = {'LexAact1', 'LexAact2', 'LexAact3', 'LexAact4'};
	elseif strcmp(varargin{1},'variableformulas'),
		% Return variable formulas in cell-array
		output = {'kon1/(kon1+LexA1)', 'kon2/(kon2+LexA2)', 'kon3/(kon3+LexA3)', 'kon4/(kon4+LexA4)'};
	else
		error('Wrong input arguments! Please read the help text to the ODE file.');
	end
	output = output(:);
	return
elseif nargin == 2,
	time = varargin{1};
	statevector = varargin{2};
elseif nargin == 3,
	time = varargin{1};
	statevector = varargin{2};
	parameterValuesNew = varargin{3};
	if length(parameterValuesNew) ~= 7,
		parameterValuesNew = [];
	end
elseif nargin == 4,
	time = varargin{1};
	statevector = varargin{2};
	parameterValuesNew = varargin{4};
else
	error('Wrong input arguments! Please read the help text to the ODE file.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LexA1 = statevector(1);
LexA2 = statevector(2);
LexA3 = statevector(3);
LexA4 = statevector(4);
Mutator1 = statevector(5);
Mutator2 = statevector(6);
Mutator3 = statevector(7);
Mutator4 = statevector(8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parameterValuesNew),
	kon1 = 10;
	kon2 = 100;
	kon3 = 0.1;
	kon4 = 1;
	N = 4;
	p = 1;
	d = 0.1;
else
	kon1 = parameterValuesNew(1);
	kon2 = parameterValuesNew(2);
	kon3 = parameterValuesNew(3);
	kon4 = parameterValuesNew(4);
	N = parameterValuesNew(5);
	p = parameterValuesNew(6);
	d = parameterValuesNew(7);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LexAact1 = kon1/(kon1+LexA1);
LexAact2 = kon2/(kon2+LexA2);
LexAact3 = kon3/(kon3+LexA3);
LexAact4 = kon4/(kon4+LexA4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LexA1_dot = p+d*LexA1;
LexA2_dot = p+d*LexA2;
LexA3_dot = p+d*LexA3;
LexA4_dot = p+d*LexA4;
Mutator1_dot = p+(1/LexAact1)-d*Mutator1;
Mutator2_dot = p+(1/LexAact2)-d*Mutator2;
Mutator3_dot = p+(1/LexAact3)-d*Mutator3;
Mutator4_dot = p+(1/LexAact4)-d*Mutator4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE ODEs
output(1) = LexA1_dot;
output(2) = LexA2_dot;
output(3) = LexA3_dot;
output(4) = LexA4_dot;
output(5) = Mutator1_dot;
output(6) = Mutator2_dot;
output(7) = Mutator3_dot;
output(8) = Mutator4_dot;
% return a column vector 
output = output(:);
return


