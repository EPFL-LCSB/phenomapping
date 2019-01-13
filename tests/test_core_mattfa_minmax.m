% Tests for FVA with the mat-tfa core module
%
% USAGE:
%
%    run('test_core_mattfa_minmax.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% run settings.m before or uncomment this part or load the 
% PhenoMappingSettings.mat

% clear
% clc
% 
% this_file_directory = mfilename('fullpath');
% run(strrep(this_file_directory,'test_core_mattfa_minmax','settings.m'))



%%%%%%%%%%%%%%%%%%%%%%%
% MINMAX WITH MAT-TFA %
%%%%%%%%%%%%%%%%%%%%%%%

% Description: Get FVA (minmax) of the model in the reference conditions  
% (normally in rich medium and without any data integrated)

% inputs
precision = 5;          % nb of decimals for rounding
tagFBA = 1;             % true for FBA, false for TFA
perObj = 0;             % percentage of optimal growth
filename = strcat(modeldescription,'_PhenoMappingMinMax');

% FVA
[minmax, description, report] = checkDirectionality(model, [], ...
    tagFBA, perObj);
indRxnBlock = find(ismember(description,'blocked'));
minmax = round(minmax, precision);

save(strcat(saving_directory,filename,'.mat'));