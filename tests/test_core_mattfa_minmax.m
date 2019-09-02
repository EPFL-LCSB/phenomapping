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
% (normally in rich medium and without any data integrated).

% inputs
precision = 5;          % nb of decimals for rounding flux values
tagFBA = 0;             % true for FBA, false for TFA
perObj = 0;             % percentage of optimal growth. Here we do not fix any growth requirement (perObj = 0). By defining perObj>0 we set a lower bound in the biomass reaction.
filename = strcat(modeldescription,'_PhenoMappingMinMax');

% FVA
[minmax, description, report] = checkDirectionality(model, [], ...
    tagFBA, perObj);

% identify blocked reactions
indRxnBlock = find(ismember(description,'blocked'));

% round flux values
minmax = round(minmax, precision);

save(strcat(saving_directory,filename,'.mat'));