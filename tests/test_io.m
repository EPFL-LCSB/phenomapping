% Tests for the phenomapping io module
%
% USAGE:
%
%    run('test_io.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%


% run settings.m before or uncomment this part or load the
% PhenoMappingSettings.mat

% clear
% clc
%
% this_file_directory = mfilename('fullpath');
% run(strrep(this_file_directory,'test_io','settings.m'))

tag_save_substrates = 0;
tag_save_substrates_joint = 0;
tag_save_secretions = 0;
tag_save_secretions_joint = 0;
tag_save_metabolomics = 0;
tag_save_transcriptomics = 1;

if ~exist('modeldescription','var') || ~exist('saving_directory','var')
    [modeldescription, saving_directory] = retrievePaths([], []);
    load(strcat(saving_directory,modeldescription,'_PhenoMappingSettings.mat'));
end

%%%%%%%%%%%%%%%%%%%%%%
% SUBSTRATE ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%

if tag_save_substrates || tag_save_secretions
    if tag_save_substrates
        filename = strcat(modeldescription,'_PhenoMappingSubstrates');
    elseif tag_save_secretions
        filename = strcat(modeldescription,'_PhenoMappingSecretions');
    end
    r1 = load(strcat(saving_directory,filename,'.mat'));
    
    % Extract info about composition of the IMMs
    [immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
        immOutput.statsClass] = extractInfoIMMDPs(r1.modelmilp, r1.DPsimm, ...
        r1.modelmilp.indUSE);
    
    immOutput.summary = [strrep(immOutput.Mets,',',''), ...
        immOutput.drainClass, num2cell(immOutput.StatsMets)];
    
    writeData(strcat(saving_directory,filename,'.csv'), immOutput.summary,...
        '%s\t%s\t%s\t%s\t%i', {'Drain ID', 'Met ID', 'Met name', ...
        'Classification based on IMM', 'Appearance in IMM'}, ...
        '%s\t%s\t%s\t%s\t%s');
    
    % Extract info: substrates linked to essentiality of the IMMs
    exportPhenoMappingInfo(r1.essIMMaddToIRM, r1.subsToGenes, 'IMM', ...
        phenotypes_directory, phenotypes_description,...
        strcat(saving_directory,filename,'_ess2sub4imm'));
end



if tag_save_substrates_joint || tag_save_secretions_joint
    if tag_save_substrates_joint
        filename = strcat(modeldescription,'_PhenoMappingSubstratesJoint');
    elseif tag_save_secretions_joint
        filename = strcat(modeldescription,'_PhenoMappingSecretionsJoint');
    end
    r1 = load(strcat(saving_directory,filename,'.mat'));
    
    % Extract info: substrates linked to essentiality of the IMMs
    exportPhenoMappingInfo(r1.essIMMaddToIRMJoint, r1.subsToGenesJoint, ...
        'IMMJoint', phenotypes_directory, phenotypes_description,...
        strcat(saving_directory,filename,'_ess2sub4immJoint'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% METABOLOMICS ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%

if tag_save_metabolomics
    filename = strcat(modeldescription,'_PhenoMappingMetabolomics');
    r2 = load(strcat(saving_directory,filename,'.mat'));
    
    % Extract info: metabolite concentrations linked to essentiality
    exportPhenoMappingInfo(r2.addEssMetab, r2.botMets, 'Metabolites',...
        phenotypes_directory, phenotypes_description, ...
        strcat(saving_directory,filename,'_metID'));
    
    exportPhenoMappingInfo(r2.addEssMetab, r2.botMetNames, ...
        'Metabolites', phenotypes_directory, phenotypes_description, ...
        strcat(saving_directory,filename,'_metNames'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSCRIPTOMICS ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tag_save_transcriptomics
    filename = strcat(modeldescription,'_PhenoMappingTranscriptomics');
    r3 = load(strcat(saving_directory,filename,'.mat'));
    
    % Extract info: metabolite concentrations linked to essentiality
    exportPhenoMappingInfo(r3.addEssGenesExp, r3.botRxnLevels, 'Rxn Levels',...
        phenotypes_directory, phenotypes_description, ...
        strcat(saving_directory,filename));
end
