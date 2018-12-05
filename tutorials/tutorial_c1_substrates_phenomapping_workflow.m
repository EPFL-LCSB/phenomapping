% Example script for PhenoMapping analysis based on substrate availability
% (uptakes: with the in silico minimal media analysis (IMM), and secretion: 
%  with the in silico minimal secretion analysis (IMS))
cd ..
%% inputs
% model should have been loaded for analysis - if it was not uncomment this
% modeldescription = 'iPbe liver';
% modelPath = '/phenomapping/models/pbe/tipbe2_liver.mat';
% load(modelPath)
% model = tipbe_liver;

% sugested inputs for imm analysis
flagUpt = 1; % true for analysis of IMM, false for analysis of IMS
minObj = 0.05; % minimal required growth
maxObj = 10; % upper bound in growth (leave unconstrained)
NumAlt = 1;%5000; % number of alternatives
drainsForiMM = {}; % apply IMM in all drains
tagMin = 1; % additional constrain to avoid generating suboptimal solutions
rxnNoThermo = []; % empty to keep analysis on structure of the model as it is
ReactionDB = load('/phenomapping/ext/matTFA/thermoDatabases/thermo_data.mat');
ReactionDB = ReactionDB.DB_AlbertyUpdate;
metabData = []; % no metabolomics data
time = []; % time limit for optimization in seconds, if empty we do not set a time limit
filename = strcat(modeldescription,'_imm');
essThr = 0.1; % essentiality threshold

%% IMM analysis

% prepare model for IMM analysis and define milp problem
[modelmilp, drains, modelpre] = analysisIMM(model, flagUpt, minObj, maxObj, ...
    drainsForiMM, rxnNoThermo, ReactionDB, metabData);

% get alternative solutions
[DPs, model4DP] = findDPMinMets(modelmilp, NumAlt, modelmilp.indUSE, ...
    time, tagMin, filename);
save(strcat('tmpresults/',filename,'_final.mat'));

% Extract info about composition of the IMMs
[immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
    immOutput.statsClass] = extractInfoIMMDPs(model4DP, DPs, model4DP.indUSE);
immOutput.summary = [strrep(immOutput.Mets,',',''), immOutput.drainClass, num2cell(immOutput.StatsMets)];
writeData(strcat('tmpresults/',filename,'.csv'),immOutput.summary,...
    '%s\t%s\t%s\t%s\t%i', {'Drain ID', 'Met ID', 'Met name', ...
    'Classification based on IMM', 'Appearance in IMM'}, '%s\t%s\t%s\t%s\t%s');


%% Get essentiality at the rich medium
grRatioGeneTFA = thermoSingleGeneDeletion(modelpre, 'TFA', modelpre.genes, 0, 0, 0, essThr, modelpre.indNF);
grRatioGeneTFA(isnan(grRatioGeneTFA)) = 0; %by default NaN is considered an essential KO
essIRMtfa = modelpre.genes(grRatioGeneTFA<essThr);

%% Option 1: get essentiality at the IMMs and link it to the substrates missing at each IMM
[essIMMfba, solOptimmfba] = getEssGeneIMM(modelpre, DPs, modelmilp, 'FBA', essThr, [], flagUpt, [], filename);
[essIMMtfa, solOptimmtfa] = getEssGeneIMM(modelpre, DPs, modelmilp, 'TFA', essThr, essIMMfba, flagUpt, [], filename);
essIMMfbatfa = getOverlapSets(essIMMfba,essIMMtfa);
save(strcat('tmpresults/',filename,'_ess_final.mat'));

[subsToGenes, essIMMaddToIRM] = linkEssGeneIMM2Subs(modelmilp, essIMMfbatfa, DPs, modelpre, essThr, essIRMtfa, NumAlt, time, [], filename);
save(strcat('tmpresults/',filename,'_ess_sub_final.mat'));

% Extract info about substrates linked to essentiality of the joint IMMs
exportIMM2ess2subInfo(essIMMaddToIRM,subsToGenes,strcat(filename,'_ess_sub_final'));


%% Option 2: get essentiality at the jointIMM (here you might add or take 
% out substrates based on you knowledge on the uptakes of your organism at 
% the life-stage of study) and link it to the substrates missing at each IMM
% extract information from IMMs
jointIMM = immOutput.Mets(immOutput.StatsMets>0.1,:);

[essJointIMMtfa, solOptjimmtfa] = getEssGeneIMM(modelpre, DPs, modelmilp, 'TFA', essThr, {essIRMtfa}, flagUpt, jointIMM(:,1), filename);
[subsToGenesJoint, essIMMaddToIRMJoint] = linkEssGeneIMM2Subs(modelmilp, essJointIMMtfa, DPs, modelpre, essThr, essIRMtfa, NumAlt, time, jointIMM(:,1), filename);
save(strcat('tmpresults/',filename,'_ess_subjoint_final.mat'));

% Extract info about substrates linked to essentiality of the joint IMMs
exportIMM2ess2subInfo(essIMMaddToIRMJoint,subsToGenesJoint,strcat(filename,'_ess_subjoint_final'));
