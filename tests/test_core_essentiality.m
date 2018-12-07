% Tests for the phenomapping essentiality core module
%
% USAGE:
%
%    run('test_core_essentiality.m')
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
% run(strrep(this_file_directory,'test_core_essentiality','settings.m'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENE ESSENTIALITY WITH TFA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description: Get essentiality of the model in the reference conditions  
% (normally in rich medium and without any data integrated)

% inputs
grRate = 0.05;          % optimal growth rate
essThr = 0.1;           % essentiality threshold
filename = strcat(modeldescription,'_PhenoMappingEssentiality');

% analysis
[~, grRateGeneTFA] = thermoSingleGeneDeletion(model, 'TFA', ...
    model.genes, 0, 0, 0, essThr, model.indNF);
grRateGeneTFA(isnan(grRateGeneTFA)) = 0; %by default NaN is considered an essential KO
essTFAref = model.genes(grRateGeneTFA < essThr*grRate);

save(strcat(saving_directory,filename,'.mat'));