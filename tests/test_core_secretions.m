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
grRate = 0.05;          % optimal growth rate
essThr = 0.1;           % essentiality threshold
flagUpt = 0;            % true for analysis of IMM, false for analysis of IMS
minObj = essThr*grRate; % minimal required growth
maxObj = 10;            % upper bound in growth (leave unconstrained)
drainsForiMM = {};      % apply IMM in all drains
metabData = [];         % no metabolomics data
NumAlt = 1;             % 1 for test! suggested 5000; % number of alternatives
time = [];              % time limit for optimization in seconds, if empty we do not set a time limit
tagMin = 1;             % additional constrain to avoid generating suboptimal solutions
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
