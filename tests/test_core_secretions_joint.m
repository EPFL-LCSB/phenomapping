% Tests for the phenomapping core substrates module
%
% USAGE:
%
%    run('test_core_secretions_joint.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% run test_core_secretions.m before or uncomment this part or load the
% PhenoMappingSecretions.mat

% clear
% clc
% 
% this_file_directory = mfilename('fullpath');
% run(strrep(this_file_directory,'_secretions_joint','_substrates.m'))



%%%%%%%%%%%%%%%%%%%%%%
% SUBSTRATE ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%

% inputs
grRate = 0.05;          % optimal growth rate
essThr = 0.1;           % essentiality threshold
flagUpt = 0;            % true for analysis of IMM, false for analysis of IMS
minObj = essThr*grRate; % minimal required growth
NumAlt = 20;            % 20 for test! suggested 5000; % number of alternatives
time = [];              % time limit for optimization in seconds, if empty we do not set a time limit

filename = strcat(modeldescription,'_PhenoMappingSecretionsJoint');

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
