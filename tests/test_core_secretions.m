% Tests for the phenomapping core substrates module
%
% USAGE:
%
%    run('test_core_secretions.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% run test_core_essentiality.m before or uncomment this part or load the
% PhenoMappingEssentiality.mat

% clear
% clc
% 
% this_file_directory = mfilename('fullpath');
% run(strrep(this_file_directory,'_secretions','_essentiality.m'))



%%%%%%%%%%%%%%%%%%%%%%
% SUBSTRATE ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%

% inputs
grRate = 0.05;          % optimal growth rate (this can be obtained optimizing for growth)
essThr = 0.1;           % essentiality threshold (% of optimal growth to be defined, see minObj). If a knockout leads to grow below this threshold the gene will be considered essential for growth.
minObj = essThr*grRate; % minimal required growth
flagUpt = 0;            % true for analysis of in silico minimal media (IMM). False for analysis of in silico minimal secretion (IMS)
maxObj = 10;            % upper bound in growth (not applied here since 10 is much higher than the growth yield of any GEM)
drainsForiMM = {};      % names of drains/exchanges used for IMM analysis. Empty means by default it will search for minimal media accross all drains/exchanges
metabData = [];         % provide metabolomics data to integrate in this analysis. Empty if none.
NumAlt = 1;             % number of alternative minimal media to find. We define 1 for you to test that the pipeline works. But it is suggested to define 5000 or more to make sure you identify all alternatives
time = [];              % time limit for optimization in seconds. If empty we do not set a time limit
tagMin = 1;             % additional constrain to avoid generating suboptimal solutions, meaning that if 1 we will not identify media that is not minimal 
filename = strcat(modeldescription,'_PhenoMappingSecretions');

% IMM analysis
% Define milp problem
modelmilp = analysisIMM(model, flagUpt, minObj, maxObj, drainsForiMM, ...
     metabData);

% Get alternative solutions
DPsimm = findDPMax(modelmilp, NumAlt, modelmilp.indUSE, ...
    time, tagMin, strcat(saving_directory,filename));

% Get essentiality at the IMMs and link it to the substrates missing at 
% each IMM
essIMMfba = getEssGeneIMM(model, DPsimm, modelmilp, 'FBA', minObj, [], ...
     flagUpt, [], strcat(saving_directory,filename));
essIMMtfa = getEssGeneIMM(model, DPsimm, modelmilp, 'TFA', minObj, ...
    essIMMfba, flagUpt, [], strcat(saving_directory,filename));
essIMMfbatfa = getOverlapSets(essIMMfba, essIMMtfa);

[subsToGenes, essIMMaddToIRM] = linkEssGeneIMM2Subs(modelmilp, ...
    essIMMfbatfa, DPsimm, model, minObj, essTFAref, NumAlt, time, [], ...
    strcat(saving_directory,filename));

save(strcat(saving_directory,filename,'.mat'));

elimIntermPMFile(saving_directory,filename)
