% Tests for the phenomapping core transcriptomics module
%
% USAGE:
%
%    run('test_core_transcriptomics.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% run test_core_mattfa_minmax.m and test_core_essentiality.m before or 
% uncomment this part or load the PhenoMappingMinMax.mat and the 
% PhenoMappingEssentiality.mat

% clear
% clc
% 
% this_file_directory = mfilename('fullpath');
% % run(strrep(this_file_directory,'_transcriptomics','_mattfa_minmax.m'))
% % run(strrep(this_file_directory,'_transcriptomics','_essentiality.m'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSCRIPTOMICS ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs
grRate = 0.05;          % optimal growth rate (this can be obtained optimizing for growth)
essThr = 0.1;           % essentiality threshold (% of optimal growth to be defined, see minObj). If a knockout leads to grow below this threshold the gene will be considered essential for growth.
minObj = essThr*grRate; % minimal required growth
NumAlt = 1;             % number of alternative minimal media to find. We define 1 for you to test that the pipeline works. But it is suggested to define 5000 or more to make sure you identify all alternatives
lowPvalue = 25;         % percentile to assign lowly expressend genes
highPvalue = 75;        % percentile to assign highly expressend genes
percent_l = 2E-3;       % determines upper bound on fluxes associated to lowly expressed genes
percent_h = 2E-5;       % determines lower bound on fluxes associated to highly expressed genes
selectAlt = [];         % expression profile that we desire to fix, empty just requires maximum consistency score (MCS)
flagPlot = 1;           % true to see plot of gene expression levels
geneKO = {};            % provide gene ids to knock out
rxnLevelOUT = model.rxns(indRxnBlock); % rxn ids that we should no consider for the TEX-FBA analyses. Here we disregard the blocked reactions based on the minmax.
indInt = [];            % indexes of UP and DOWN rxns (once TEX-FBA constraints are integrated) among which we want to identify bottleneck reactions. Empty means we will go with the default option: looking among all TEX-FBA reaction constraints.
filename = strcat(modeldescription,'_PhenoMappingTranscriptomics');


% TEX-FBA analysis
% load levelGenes.mat
load(transcriptomics_directory);

% Integrate gene expression constraints
[expModelAlt, solution, newgrRules, expModel] = integrateGeneExp ...
    (model, levelGenes, minObj, lowPvalue, highPvalue, percent_l, ...
    percent_h, NumAlt, selectAlt, flagPlot, minmax, geneKO, rxnLevelOUT);

% Get essentiality at each expression profile
[addEssGenesExp, essGenesExp, grRateExp] = getEssGeneExp(expModel, ...
    solution, model, minObj, essTFAref);

% Link gene essentiality at each expression profile to rxn levels
% responsible for it
botRxnLevels = linkEssGene2Exp(expModel, addEssGenesExp, solution, ...
    indInt, minObj, NumAlt, model.indNF, strcat(saving_directory,...
    filename,'_DPs'));

save(strcat(saving_directory,filename,'.mat'));

elimIntermPMFile(saving_directory,filename)