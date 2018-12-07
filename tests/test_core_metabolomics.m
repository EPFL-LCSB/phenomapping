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
grRate = 0.05;          % optimal growth rate
essThr = 0.1;           % essentiality threshold
minObj = essThr*grRate; % minimal required growth
NumAlt = 500;           % number of alternatives
time = [];              % time limit for optimization in seconds, if empty we do not set a time limit
tagMax = 0;             % additional constrain to avoid generating suboptimal solutions
filename = strcat(modeldescription,'_PhenoMappingMetabolomics');

% Integrate metabolomics data
tmp = load(metabolomics_directory);
[modelMetab, LCcons] = prepMetabCons(model, tmp.metNames, tmp.dataLC(:,1), ...
        tmp.dataLC(:,2));
clear tmp

% Get essentiality with metabolomics data
[~, grRateGeneTFAmetab] = thermoSingleGeneDeletion(modelMetab, 'TFA',...
    modelMetab.genes, 0, 0, 0, essThr, modelMetab.indNF);
grRateGeneTFAmetab(isnan(grRateGeneTFAmetab)) = 0; %by default NaN is considered an essential KO
essTFAmetab = modelMetab.genes(grRateGeneTFAmetab < minObj);

addEssMetab = essTFAmetab(~ismember(essTFAmetab,essTFAref));

% Get bottleneck metabolites
bottleneckMets = cell(length(addEssMetab),1);
bottleneckMetNames = cell(length(addEssMetab),1);

for i = 1:length(addEssMetab)
    modelDel = thermoDeleteModelGenes(model, addEssMetab{i});
    [bottleneckMets{i}, ~, bottleneckMetNames{i}] = getBotNeckMets(...
        modelDel, LCcons, minObj, NumAlt, time, tagMax, ...
        strcat(saving_directory,filename));
end

save(strcat(saving_directory,filename,'.mat'));

elimIntermPMFile(saving_directory,filename)
