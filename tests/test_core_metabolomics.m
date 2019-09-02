% Tests for the phenomapping core metabolomics module
%
% USAGE:
%
%    run('test_core_metabolomics.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% run test_core_essentiality.m before or uncomment this part  or load the
% PhenoMappingEssentiality.mat

% clear
% clc
%
% this_file_directory = mfilename('fullpath');
% % run(strrep(this_file_directory,'_metabolomics','_essentiality.m'))



%%%%%%%%%%%%%%%%%%%%%%%%%
% METABOLOMICS ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs
grRate = 0.05;          % optimal growth rate (this can be obtained optimizing for growth)
essThr = 0.1;           % essentiality threshold (% of optimal growth to be defined, see minObj). If a knockout leads to grow below this threshold the gene will be considered essential for growth.
minObj = essThr*grRate; % minimal required growth
NumAlt = 500;           % number of alternative minimal media to find. It is suggested to define a big number like 500.
tagMax = 0;             % additional constrain to avoid generating suboptimal solutions, meaning that if 1 we will not identify metabolite sets that is not maximal.
filename = strcat(modeldescription,'_PhenoMappingMetabolomics');

% Integrate metabolomics data
tmp = load(metabolomics_directory);
[modelMetab, LCcons, tagFeas] = prepMetabCons(model, tmp.metNames, ...
    tmp.dataLC(:,1), tmp.dataLC(:,2));
clear tmp

% Get essentiality with metabolomics data
[~, grRateGeneTFAmetab] = thermoSingleGeneDeletion(modelMetab, ...
    'TFA', modelMetab.genes, 0, 0, 0, essThr, modelMetab.indNF);
grRateGeneTFAmetab(isnan(grRateGeneTFAmetab)) = 0;
essTFAmetab = modelMetab.genes(grRateGeneTFAmetab < minObj);

% Get bottleneck metabolites
[botMets, botMetNames, addEssMetab] = ...
    linkEssGeneMetab2Mets(model, LCcons, minObj, essTFAref, ...
    essTFAmetab, NumAlt, time, tagMax, strcat(saving_directory,filename));

save(strcat(saving_directory,filename,'.mat'));

elimIntermPMFile(saving_directory,filename)
