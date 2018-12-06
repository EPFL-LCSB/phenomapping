% Example script for PhenoMapping analysis based on substrate availability
% (uptakes: with the in silico minimal media analysis (IMM), and secretion: 
%  with the in silico minimal secretion analysis (IMS))
cd(pathToPhenoMapping)
%% Inputs
% the model should have been loaded for this analysis and gone through the 
% initTestPhenoMappingModel test (as described in the 
% master_tutorial_phenomapping_workflow.m)

% sugested inputs for PhenoMapping analysis of substrates
flagUpt = 1;            % true for analysis of IMM, false for analysis of IMS
minObj = 0.05;          % minimal required growth
maxObj = 10;            % upper bound in growth (leave unconstrained)
drainsForiMM = {};      % apply IMM in all drains
metabData = [];         % no metabolomics data
NumAlt = 1;             % 1 for test! suggested 5000; % number of alternatives
time = [];              % time limit for optimization in seconds, if empty we do not set a time limit
tagMin = 1;             % additional constrain to avoid generating suboptimal solutions
essThr = 0.1;           % essentiality threshold

filename = strcat(modeldescription,'_PhenoMappingSubstrates');
if ~exist('pathToPhenoData','var')
    pathToPhenoData = [];
end
if ~exist('phenDesc','var')
    phenDesc = pathToPhenoData;
end
if ~exist('pathToSave','var')
    pathToSave = 'tmpresults/';
end

%% IMM analysis
% Define milp problem
[modelmilp, drains] = analysisIMM(model, flagUpt, minObj, maxObj, ...
    drainsForiMM, metabData);

% Get alternative solutions
[DPs, model4DP] = findDPMax(modelmilp, NumAlt, modelmilp.indUSE, ...
    time, tagMin, strcat(pathToSave,filename));

save(strcat(pathToSave,filename,'.mat'));

% Extract info about composition of the IMMs
[immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
    immOutput.statsClass] = extractInfoIMMDPs(model4DP, DPs, ...
    model4DP.indUSE);
immOutput.summary = [strrep(immOutput.Mets,',',''), ...
    immOutput.drainClass, num2cell(immOutput.StatsMets)];
writeData(strcat(pathToSave,filename,'.csv'), immOutput.summary,...
    '%s\t%s\t%s\t%s\t%i', {'Drain ID', 'Met ID', 'Met name', ...
    'Classification based on IMM', 'Appearance in IMM'}, ...
    '%s\t%s\t%s\t%s\t%s');


%% Get essentiality at the rich medium
if ~exist('essTFAref','var') || isempty(essTFAref)
    [~, grRateGeneTFA] = thermoSingleGeneDeletion(model, 'TFA', ...
        model.genes, 0, 0, 0, essThr, model.indNF);
    grRateGeneTFA(isnan(grRateGeneTFA)) = 0; %by default NaN is considered an essential KO
    essTFAref = model.genes(grRateGeneTFA < minObj);
end

%% Option 1: get essentiality at the IMMs and link it to the substrates missing at each IMM

[essIMMfba, solOptimmfba] = getEssGeneIMM(model, DPs, modelmilp, 'FBA', ...
    essThr, [], flagUpt, [], strcat(pathToSave,filename));
[essIMMtfa, solOptimmtfa] = getEssGeneIMM(model, DPs, modelmilp, 'TFA', ...
    essThr, essIMMfba, flagUpt, [], strcat(pathToSave,filename));
essIMMfbatfa = getOverlapSets(essIMMfba, essIMMtfa);

save(strcat(pathToSave,filename,'.mat'));

[subsToGenes, essIMMaddToIRM] = linkEssGeneIMM2Subs(modelmilp, ...
    essIMMfbatfa, DPs, model, essThr, essIRMtfa, NumAlt, time, [], ...
    strcat(pathToSave,filename));

save(strcat(pathToSave,filename,'.mat'));

% Extract info: substrates linked to essentiality of the IMMs
exportPhenoMappingInfo(essIMMaddToIRM, subsToGenes, 'IMM', ...
    pathToPhenoData, phenDesc,strcat(pathToSave,filename,'_ess2sub4imm'));

%% Option 2: get essentiality at the jointIMM (here you might add or take 
% out substrates based on you knowledge on the uptakes of your organism at 
% the life-stage of study) and link it to the substrates missing in the 
% joinIMM
jointIMM = immOutput.Mets(immOutput.StatsMets > 0.1, :);

[essJointIMMtfa, solOptjimmtfa] = getEssGeneIMM(model, DPs, modelmilp, ...
    'TFA', essThr, {essIRMtfa}, flagUpt, jointIMM(:,1), ...
    strcat(pathToSave,filename));
[subsToGenesJoint, essIMMaddToIRMJoint] = linkEssGeneIMM2Subs(...
    modelmilp, essJointIMMtfa, DPs, model, essThr, essIRMtfa, NumAlt, ...
    time, jointIMM(:,1), strcat(pathToSave,filename));

save(strcat(pathToSave,filename,'.mat'));

% Extract info: substrates linked to essentiality of the joint IMMs
exportPhenoMappingInfo(essIMMaddToIRMJoint, subsToGenesJoint, ...
    'jointIMM', pathToPhenoData, phenDesc, strcat(pathToSave, ...
    filename, '_ess2sub4jointimm'));
