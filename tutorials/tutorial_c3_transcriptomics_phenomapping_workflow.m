% Example script for PhenoMapping analysis based on absolute gene
% expression data
cd(pathToPhenoMapping)
%% inputs
% the model should have been loaded for this analysis and gone through the
% initTestPhenoMappingModel test (as described in the
% master_tutorial_phenomapping_workflow.m)

% sugested inputs for PhenoMapping analysis of gene expression/reaction
% levels
minObj = 0.05;          % minimal required growth
NumAlt = 500;           % number of alternatives
essThr = 0.1;           % essentiality threshold
selectAlt = 0;          % zero to define lb on CS and not to fix the profile of an alternative
lowPvalue = 20;         % based on the distribution: genes that are considered to be lowly expressed (see TEXFBA paper for details)
highPvalue = 91.8;      % based on the distribution: genes that are considered to be highly expressed (see TEXFBA paper for details)
percent_h = 2E-3;       % lb flux constraint on highly expressed reaction wrt its upper bound from minmax (see TEXFBA paper for details)
percent_l = 2E-5;       % ub flux constraint on lowly expressed reaction wrt its lower bound from minmax (see TEXFBA paper for details)
flagPlot = 1;           % True to plot the distribution of lowly and highly expressed genes based on parameters provided
rxnLevelOUT = {};       % List of reactions for which the levels (derived from the gene expression levels) should not be kept
geneKO = {};            % List of genes for KO
minmax = runMinMax;     % Run minmax

filename = strcat(modeldescription,'_PhenoMappingTranscriptomics');
if ~exist('TransPath','var')
    TransPath = [];
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
if isempty(TransPath)
    while isempty(TransPath)
        TransPath = input('Please provide the path to the metabolomics data (metNames+dataLC or LCcons) and press enter\n... ','s');
    end
end
tmp = load(TransPath);
if isfield(tmp,'levelGenes')
    levelGenes = tmp.levelGenes;
else
    error('the transcriptomics data should be saved in a numeric array called "levelGenes" (see generateGeneLevels for details)')
end

%% Integration of transcriptomics data into the model with TEXFBA

[expModel, solution, newgrRules, expModelm] = integrateGeneExp ...
    (model, levelGenes, minObj, lowPvalue, highPvalue, percent_l, ...
    percent_h, NumAlt, selectAlt, flagPlot, minmax, geneKO, rxnLevelOUT);

save(strcat(pathToSave,filename,'.mat'));

%% Option 1: Get essentiality with transcriptomics data considering there 
% can be regulation of expression between isoenzymes

% for MCS 
[~, grRateGeneTFAge] = thermoSingleGeneDeletion(expModel, 'TFA', ...
    expModel.genes, 0, 0, 0, essThr, expModel.indNF);
grRateGeneTFAge(isnan(grRateGeneTFAge)) = 0; %by default NaN is considered an essential KO
essTFAge = expModel.genes(grRateGeneTFAge < minObj);

if ~exist('essTFAref','var') || isempty(essTFAref)
    [~, grRateGeneTFA] = thermoSingleGeneDeletion(model, 'TFA', ...
        model.genes, 0, 0, 0, essThr, model.indNF);
    grRateGeneTFA(isnan(grRateGeneTFA)) = 0; %by default NaN is considered an essential KO
    essTFAref = model.genes(grRateGeneTFA < minObj);
end
addEssMetab = essTFAge(~ismember(essTFAge,essTFAref));

save(strcat(pathToSave,filename,'.mat'));

% for each transcriptomic profile
indUPDOWN = [expModelm.ExpInfo.indUP; expModelm.ExpInfo.indDOWN];
essGenesTransAlt = cell(size(solution.z_matrix,2),1);
grRateGeneTFAgeAlt = cell(size(solution.z_matrix,2),1);

for sAi = 1:size(solution.z_matrix,2) % select alternative to integrate
    expModelAlt = integrateConsProfile(expModelm, solution, indUPDOWN, ...
        sAi);
    [~, grRateGeneTFAgeAlt{sAi}] = thermoSingleGeneDeletion(expModelAlt,...
        'TFA', expModelAlt.genes, 0, 0, 0, essThr, expModelAlt.indNF);
    grRateGeneTFAgeAlt{sAi}(isnan(grRateGeneTFAgeAlt{sAi})) = 0;
    essGenesTransAlt{sAi,1} = expModelm.genes(grRateGeneTFAgeAlt{sAi} < ...
        minObj);
end

save(strcat(pathToSave,filename,'.mat'));

%% Get bottleneck reaction levels

allToRelax = cell(1,size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2));
addEssGenesAs = cell(1,size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2));

for sAi = 1:size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2) % select alternative to integrate
    newess = mm.essGenes_expModelAlt{sAi,col};
    addEssGenes = newess(not(ismember(newess,essorig))); %additional genes identified with gene expression data integrated vs previous version (TFA tipbe2 liver)
    addEssGenesAs{1,sAi} = addEssGenes;    
    indObj = indUPDOWN(solution.sol_matrix(indUPDOWN,sAi)>0.9);
    for i = 1:length(addEssGenes)
        modeli = thermoDeleteModelGenes(expModelm,addEssGenes(i));
        [genesToKeep,genesToRelax] = getBotGenes(modeli,indObj,grRate,numAlt,indNF);
        allToRelax{i,sAi} = genesToRelax;
    end
end


% Extract info: reactions levels linked to essentiality
exportPhenoMappingInfo(addEssMetab,bottleneckMets,'Metabolites',pathToPhenoData,phenDesc,strcat(filename,'_metID'));
exportPhenoMappingInfo(addEssMetab,bottleneckMetNames,'Metabolites',pathToPhenoData,phenDesc,strcat(filename,'_metNames'));

save(strcat(pathToSave,filename,'.mat'));

%% Option 2: Get essentiality with transcriptomics data considering there
% cannot be regulation of expression between isoenzymes (lowly expressed
% genes are knocked out)

% for MCS 
[~, essTFAgeNOREG] = thermoSingleGeneDeletionGeneExpNewGPR(expModel, ...
    minObjCHECK, essThr, 0);
addEssMetabNOREG = essTFAgeNOREG(~ismember(essTFAgeNOREG,essTFAref));

save(strcat(pathToSave,filename,'.mat'));

% for each transcriptomic profile
essGenesTransAltNOREG = cell(size(solution.z_matrix,2),1);

for sAi = 1:size(solution.z_matrix,2) % select alternative to integrate
    expModelAlt = integrateConsProfile(expModelm, solution, indUPDOWN, ...
        sAi);
    [~, essGenesTransAltNOREG{sAi,1}, ~] = ...
        thermoSingleGeneDeletionGeneExpNewGPR(expModelAlt, minObjCHECK, ...
        essThr, 0);
end

save(strcat(pathToSave,filename,'.mat'));