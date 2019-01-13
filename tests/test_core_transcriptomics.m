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
indTransDir = 1;        % index in transcriptomics_directory
grRate = 0.05;          % optimal growth rate
essThr = 0.1;           % essentiality threshold
minObj = essThr*grRate; % minimal required growth
lowPvalue = 25;
highPvalue = 75;
percent_l = 2E-3;
percent_h = 2E-5;
NumAlt = 1;             % 1 for test! suggested 5000; % number of alternatives
selectAlt = [];         % expression profile that we desire to fix, empty just requires MCS
flagPlot = 1;           % true to see plot of gene expression levels
geneKO = {};            % 
rxnLevelOUT = model.rxns(indRxnBlock);
indInt = [];
filename = strcat(modeldescription,'_PhenoMappingTranscriptomics');


% TEX-FBA analysis
% load levelGenes.mat
load(transcriptomics_directory{indTransDir});

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
    indInt, minObj, NumAlt, model.indNF, filename);

save(strcat(saving_directory,filename,'.mat'));

elimIntermPMFile(saving_directory,filename)