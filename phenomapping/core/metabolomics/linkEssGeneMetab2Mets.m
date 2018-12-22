function [botMets, botMetNames, addEssMetab] = ...
    linkEssGeneMetab2Mets(model, LCcons, minObj, essTFAref, ...
    essTFAmetab, NumAlt, time, tagMax, filename)
% Link essentiality of a gene to constraints from metabolomics data within
% TFA
%
% USAGE:
%
%       [botMets, botMetNames, addEssMetab] = linkEssGeneMetab2Mets(model, LCcons, minObj, essTFAref, essTFAmetab, NumAlt, time, tagMax, filename)
%
% INPUTS:
%    model:           FEASIBLE model with TFA structure
%    LCcons:          Default structure to integrate concentrations in TFA
%
% OPTIONAL INPUTS:
%    minObj:          min growth to achieve (default = 0.05)
%    essTFAref:       List of genes essential without metabolomics data
%                     integrated - reference (default = calculate)
%    essTFAmetab      List of genes essential with metabolomics data
%                     integrated - new condition (default = calculate)
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    time:            Time in sec after which simulation will be stopped
%                     (default = empty / no time provided)
%    tagMax:          True to generate alternative solution of minimal size
%                     only. False to generate up to the number of
%                     alternatives provided in NumAlt (default = false)
%    filename:        Name used to save DPs (default = 'PhenoMappingDPMin')
%
% OUTPUTS:
%    botMets:         Cell array with metabolites (met IDs) found in each 
%                     alternative
%    botMetNames:     Cell array with metabolites (met names) found in each 
%                     alternative
%    addEssMetab:     Essential genes with metabolomics on the top of 
%                     essential genes without metabolomics (reference
%                     condition)
%
% .. Author:
% Anush Chiappino-Pepe 2016
% 

if (nargin < 3)
    minObj = 0.05;
end
if (nargin < 4)
    [~, grRateGeneTFA] = thermoSingleGeneDeletion(model, 'TFA', ...
        model.genes, 0, 0, 0, 0, model.indNF);
    grRateGeneTFA(isnan(grRateGeneTFA)) = 0; %by default NaN is considered an essential KO
    essTFAref = model.genes(grRateGeneTFA < minObj);
end
if (nargin < 5)
    modelMetab = loadConstraints(model, LCcons);
    sol = optimizeThermoModel(modelMetab);
    if isnan(sol.val) || isempty(sol.val) || (sol.val<1E-3)
        error('You need to generate a reduced metabolomics data set for integration. Check tutorial_issues.m');
    end
    [~, grRateGeneTFAmetab] = thermoSingleGeneDeletion(modelMetab, 'TFA',...
        modelMetab.genes, 0, 0, 0, 0, modelMetab.indNF);
    grRateGeneTFAmetab(isnan(grRateGeneTFAmetab)) = 0; %by default NaN is considered an essential KO
    essTFAmetab = modelMetab.genes(grRateGeneTFAmetab < minObj);
end
if (nargin < 6)
    NumAlt = 1;
end
if (nargin < 7)
    time = [];
end
if (nargin < 8)
    tagMax = 0;
end
if (nargin < 9)
    filename = 'PhenoMappingDPMin';
end


addEssMetab = essTFAmetab(~ismember(essTFAmetab,essTFAref));

% Get bottleneck metabolites
botMets = cell(length(addEssMetab),1);
botMetNames = cell(length(addEssMetab),1);

for i = 1:length(addEssMetab)
    modelDel = thermoDeleteModelGenes(model, addEssMetab{i});
    [botMets{i}, ~, botMetNames{i}] = getBotNeckMets(...
        modelDel, LCcons, minObj, NumAlt, time, tagMax, ...
        filename);
end