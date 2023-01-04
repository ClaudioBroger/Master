function data_TetR = readdataTetR

% used to generate the one .mat file per experiment (strains measured on the same day) from the excel or facs data.

data_TetR.time = ..; % timepoint(s) of measurement


[num,txt,raw] = xlsread('data\.._.xlsx',..sheetnumber); % Read in xlsx file

data_TetR.tdh3 = ... 
data_TetR.empty = ...
data_TetR.means = ...; 
data_TetR.std = ...;

modulenames{1} = 'ptetCitrine-ptetTetR'; %Description name
%modulenames{2} = 'ptetCitrine-ptetTetR';

dose = ...
timepoint = ...

% .. means can for example be a matrix where each row is the DR of a module

save('data/data_TetR.mat','data_TetR','modulenames','dose','timepoint');


end