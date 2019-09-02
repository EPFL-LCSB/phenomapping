% Tests for the phenomapping core substrates module
%
% USAGE:
%
%    run('test_core_substrates_joint.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% run test_core_substrates.m before or uncomment this part or load the
% PhenoMappingSubstrates.mat

% clear
% clc
% 
% this_file_directory = mfilename('fullpath');
% run(strrep(this_file_directory,'_substrates_joint','_substrates.m'))



%%%%%%%%%%%%%%%%%%%%%%
% SUBSTRATE ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%

% inputs
grRate = 0.05;          % optimal growth rate (this can be obtained optimizing for growth)
essThr = 0.1;           % essentiality threshold (% of optimal growth to be defined, see minObj). If a knockout leads to grow below this threshold the gene will be considered essential for growth.
minObj = essThr*grRate; % minimal required growth
flagUpt = 1;            % true for analysis of in silico minimal media (IMM). False for analysis of in silico minimal secretion (IMS)
NumAlt = 20;            % number of alternative minimal media to find. We define 20 for you to test that the pipeline works. But it is suggested to define 5000 or more to make sure you identify all alternatives
time = [];              % time limit for optimization in seconds. If empty we do not set a time limit
filename = strcat(modeldescription,'_PhenoMappingSubstratesJoint');

% Get essentiality at the joint IMM and link it to the substrates missing 
% at the joint IMM
[immOutput.StatsMets, immOutput.Mets] = extractInfoIMMDPs(modelmilp, ...
    DPsimm, modelmilp.indUSE);
jointIMM = immOutput.Mets(immOutput.StatsMets > 0.1, :);

essIMMtfaJoint = getEssGeneIMM(model, DPsimm, modelmilp, 'TFA', minObj, ...
    {essTFAref}, flagUpt, jointIMM(:,1), ...
    strcat(saving_directory,filename));

[subsToGenesJoint, essIMMaddToIRMJoint] = linkEssGeneIMM2Subs(modelmilp, ...
    essIMMtfaJoint, DPsimm, model, minObj, essTFAref, NumAlt, time, ...
    jointIMM(:,1), strcat(saving_directory,filename));

save(strcat(saving_directory,filename,'.mat'));

elimIntermPMFile(saving_directory,filename)
