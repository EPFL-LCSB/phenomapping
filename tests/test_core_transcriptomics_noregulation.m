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

% run test_core_transcriptomics.m before or uncomment this part or load 
% the PhenoMappingTranscriptomics.mat

% clear
% clc
% 
% this_file_directory = mfilename('fullpath');
% % run(strrep(this_file_directory,'_noregulation',''))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSCRIPTOMICS without regulation ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs
grRate = 0.05;          % optimal growth rate
essThr = 0.1;           % essentiality threshold
minObj = essThr*grRate; % minimal required growth
filename = strcat(modeldescription,'_PhenoMappingTranscriptomicsRegulation');


% TEX-FBA analysis

% Get essentiality at each expression profile without regulation
[addEssGenesReg, essGenesReg] = getEssReg(expModel, ...
    solution, model, minObj, essTFAref);

save(strcat(saving_directory,filename,'.mat'));

elimIntermPMFile(saving_directory,filename)