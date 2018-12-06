%% Setting up the paths, solver and preparing the model for PhenoMapping
clear
clc
close all

% adding phenomapping and dependent repositories to the path
cplexPath = '/Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271'; % provide path to cplex
pathToPhenoMapping = '/Users/Anush/GIT_Folders/phenomapping'; %provide path to phenomapping repository
pathTomatTFA = '/Users/Anush/GIT_Folders/matTFA'; % provide path to matTFA repository
pathToTEXFBA = '/Users/Anush/GIT_Folders/texfba'; % provide path to TEXFBA repository

addpath(genpath(pathToPhenoMapping));
addpath(genpath(pathTomatTFA));
addpath(genpath(pathToTEXFBA));
cd(pathToPhenoMapping) % work from first folder of phenomapping as the starting path

% providing information for TFA studies
tagThermo = 1;
ReactionDBpath = '/Users/Anush/GIT_Folders/matTFA/thermoDatabases/thermo_data.mat';

% providing the model
modeldescription = 'iPbe';
modelPath = 'models/pbe/tipbe2_liver.mat';
load(modelPath)
model = tipbe_liver;

% providing path to phenotypes and omics data (if available)
pathToPhenoData = {'/Users/Anush/GIT_Folders/phenomapping/data/pbe/pbe_phenotypes_blood_Nov16.mat';
    '/Users/Anush/GIT_Folders/phenomapping/data/pbe/pbe_phenotypes_liver_Jul18.mat'};
phenDesc = {'blood pbe';'liver pbe'};
MetabPath = '/Users/Anush/GIT_Folders/phenomapping/data/pbe/allmetab_pfa_blood.mat';
TransPath = '/Users/Anush/GIT_Folders/phenomapping/data/pbe/levelGenes_pbe_liver_HepG2_mean_48h.mat';

% path to save results
pathToSave = 'tmpresults/';

% preparing model for phenomapping
[model, checkList, tagReady] = initTestPhenoMappingModel(model,cplexPath,ReactionDBpath);

save(strcat(pathToSave,'inputPhenoMapping.mat'));

%% PhenoMapping analysis for organism specific information in the GEM
% Step 1: Map essentiality to enzymatic irreversibilities defined as adhoc 
% in the model (if applicable):
% cd('tutorials')
% run('tutorial_b1_adhocirrev_phenomapping_workflow.m')


% Step 2: Map essentiality to limitation in intracellular transportability
% (if eukaryotic):
% cd('tutorials')
% run('tutorial_b2_inttransport_phenomapping_workflow.m')


% Step 3: Map essentiality to localization of enzymes (if eukaryotic):
% cd('tutorials')
% run('tutorial_b3_enyloc_phenomapping_workflow.m')


% Step 4: Map essentiality to biochemistry annotated to genome:
% cd('tutorials')
% run('tutorial_b4_biochem_phenomapping_workflow.m')

%% PhenoMapping analysis for condition specific information in the GEM
% Step 1: Map essentiality to substrate availability following the 
% PhenoMapping analysis for the medium:
cd('tutorials')
run('tutorial_c1_substrates_phenomapping_workflow.m')

% Step 2: Map essentiality to metabolite concentrations following the 
% PhenoMapping analysis for metabolomics (if available):
cd('tutorials')
run('tutorial_c2_metabolomics_phenomapping_workflow.m')

% Step 3: Map essentiality to gene expression levels following the 
% PhenoMapping analysis for transcriptomics (if available) considering or 
% not regulation of gene expression between isoenzymes:
cd('tutorials')
run('tutorial_c3_transcriptomics_phenomapping_workflow.m')

