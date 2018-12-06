% Example script for PhenoMapping analysis based on metabolite
% concentration
cd(pathToPhenoMapping)
%% Inputs
% the model should have been loaded for this analysis and gone through the
% initTestPhenoMappingModel test (as described in the
% master_tutorial_phenomapping_workflow.m)

% sugested inputs for PhenoMapping analysis of metabolite concentrations
minObj = 0.05;          % minimal required growth
NumAlt = 500;           % number of alternatives
essThr = 0.1;           % essentiality threshold
time = [];              % time limit for optimization in seconds, if empty we do not set a time limit
tagMax = 0;             % additional constrain to avoid generating suboptimal solutions

filename = strcat(modeldescription,'_PhenoMappingMetabolomics');
if ~exist('MetabPath','var')
    MetabPath = [];
end
if ~exist('pathToPhenoData','var')
    pathToPhenoData = [];
end
if ~exist('phenDesc','var')
    phenDesc = pathToPhenoData;
end
if ~exist('pathToSave','var')
    pathToSave = 'tmpresults/';
end

%% Load data
if isempty(MetabPath)
    while isempty(MetabPath)
        MetabPath = input('Please provide the path to the metabolomics data (metNames+dataLC or LCcons) and press enter\n... ','s');
    end
end
tmp = load(MetabPath);
if isfield(tmp,'metNames') && isfield(tmp,'dataLC')
    metNames = tmp.metNames;
    dataLC = tmp.dataLC;
elseif isfield(tmp,'LCcons')
    LCcons = tmp.LCcons;
else
    error('the metabolomics data should be saved in a 3-column cell called "LCcons" or in two structures "metNames" and "dataLC" (see prepMetabCons for details)')
end

%% Integration of metabolomics data into the model with MATTFA
if exist('LCcons','var')
    modelMetab = loadConstraints(model, LCcons);
elseif (exist('metNames','var') && exist('dataLC','var'))
    [modelMetab, LCcons] = prepMetabCons(model, metNames, dataLC(:,1), ...
        dataLC(:,2));
end

sol = optimizeThermoModel(modelMetab);
if isempty(sol.val) || isnan(sol.val) || sol.val<minObj
    warning('the available metabolomics data make your model infeasible. A bottleneck metabolite was selected and removed from the metabolomic set. It is suggested to manually check this issue in more detail. Please go to tutorial_issues.m')
    bottleneckMetsBasis = getBotNeckMets(model, LCcons, minObj, 1, time);
    LCconsinit = LCcons;
    LCcons = LCcons(~ismember(LCcons(:,1),strcat('LC_',...
        bottleneckMetsBasis{1,1})),:);
    modelMetab = loadConstraints(model, LCcons);
    sol = optimizeThermoModel(modelMetab);
    if isempty(sol.val) || isnan(sol.val) || sol.val<minObj
        error('Major issue with integration of metabolomics data. Please check this step carefully') % this should not happen...
    end
end

%% Get essentiality with metabolomics data
[~, grRateGeneTFAmetab] = thermoSingleGeneDeletion(modelMetab, 'TFA',...
    modelMetab.genes, 0, 0, 0, essThr, modelMetab.indNF);
grRateGeneTFAmetab(isnan(grRateGeneTFAmetab)) = 0; %by default NaN is considered an essential KO
essTFAmetab = modelMetab.genes(grRateGeneTFAmetab < minObj);

if ~exist('essTFAref','var') || isempty(essTFAref)
    [~, grRateGeneTFA] = thermoSingleGeneDeletion(model, 'TFA',...
        model.genes, 0, 0, 0, essThr, model.indNF);
    grRateGeneTFA(isnan(grRateGeneTFA)) = 0; %by default NaN is considered an essential KO
    essTFAref = model.genes(grRateGeneTFA < minObj);
end

addEssMetab = essTFAmetab(~ismember(essTFAmetab,essTFAref));

%% Get bottleneck metabolites
modelDel = cell(length(addEssMetab),1);
hasEffect = cell(length(addEssMetab),1);
constrRxnNames = cell(length(addEssMetab),1);
bottleneckMets = cell(length(addEssMetab),1);
bottleneckMetNames = cell(length(addEssMetab),1);

% Define milp problem and get bottleneck metabolites
for i = 1:length(addEssMetab)
    [modelDel{i}, hasEffect{i}, constrRxnNames{i}] = ...
        thermoDeleteModelGenes(model, addEssMetab{i});
    [bottleneckMets{i}, ~, bottleneckMetNames{i}] = getBotNeckMets(...
        modelDel{i}, LCcons, minObj, NumAlt, time, tagMax, ...
        strcat(pathToSave,filename));
end

% Extract info: metabolite concentrations linked to essentiality
exportPhenoMappingInfo(addEssMetab, bottleneckMets, 'Metabolites',...
    pathToPhenoData, phenDesc, strcat(pathToSave,filename,'_metID'));
exportPhenoMappingInfo(addEssMetab, bottleneckMetNames, 'Metabolites',...
    pathToPhenoData, phenDesc, strcat(pathToSave,filename,'_metNames'));

save(strcat(pathToSave,filename,'.mat'));