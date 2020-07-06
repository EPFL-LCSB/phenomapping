%% Initial info and contact

% This script is supposed to run in a go (after installing CPLEX, matTFA 
% and texFBA; and working from the parent phenomapping directory).

% If this is not the case: Please, contact 
% Prof. Vassily Hatzimanikatis (vassily.hatzimanikatis@epfl.ch) or
% Anush Chiappino-Pepe (anush_chiappinopepe@hms.harvard.edu)


%% Prior to working with PhenoMapping

% First, install CPLEX from: https://www.ibm.com/analytics/cplex-optimizer

% Then, install matTFA.
% Go to the same directory where phenomapping is (the parent phenomapping 
% folder that contains all subfolders: ext, models, phenomapping, tests, 
% tutorials)
% Clone the matTFA repository from GitHub by typing in your terminal:
% git clone https://github.com/EPFL-LCSB/matTFA

% Then, install texFBA.
% Go to the same directory where phenomapping is (the parent phenomapping 
% folder that contains all subfolders: ext, models, phenomapping, tests, 
% tutorials)
% Clone the texFBA repository from GitHub by typing in your terminal:
% git clone https://github.com/EPFL-LCSB/texfba


%% Setting up the paths, solver and preparing the model for PhenoMapping

% It is recommended that you remove any other file from the matlab path 
% prior to starting phenomapping

% Go inside the parent phenomapping folder and work from there.
% The parent phenomapping folder is the main repository folder that 
% contains all subfolders: ext, models, phenomapping, tests, tutorials)
% Do not run this tutorial script from other subfolders.

clear
clc
close all

% Files to set up the paths and prepare the model for PhenoMapping have
% been provided in the phenomapping/tests folder:

% Notes valid for all options:

% NOTE 1: matTFA and TEX-FBA are in the same directory as phenomapping 
% (not inside phenomapping!).
% This will be assumed by default in the "initPhenoMappingPaths.m".
% If this is not the case you will have to provide the path to matTFA and 
% TEX-FBA when you are asked to do so by "initPhenoMappingPaths.m".

% NOTE 2: You will need to paste the path to cplex when the following 
% message appears:
% "Please provide your cplex path and press enter"
% You will paste the corresponding path as follows (without " "):
% ... /Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION A: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the required repositories and solver directories to the 
% matlab path in an automatic fashion, meaning that phenomapping will 
% automatically search for matTFA and TEX-FBA in the same directory where
% the parent folder is. If these are not there you will be asked to 
% provide the paths. 
% For this option run the "settings.m" script as follows:
% run('tests/settings.m')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION B: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the required repositories and solver directories to the 
% matlab path (same as A) AND prepare the model for phenomapping analyses 
% - also done in an automatic fashion. 
% For this option run the "settings_nameofmodel.m" script as follows:
run('tests/settings_ipbeliver.m')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION C: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add manually the required repositories and solver paths to the 
% matlab path AND prepare the model for phenomapping analyses. For this
% option you can follow the backbone of the code presented below. 
% Please, note that since option C is manual, you should adapt the paths
% You should work from the parent phenomapping folder to run 
% "initTestPhenoMappingModel.m"

% Provide path to folder where your results will be saved - please, keep 
% the name "tmpresults" for this folder.
% initPhenoMappingPaths will create the "tmpresults" folder at the same 
% level as ext, models, phenomapping, tests, and tutorials subfolders.
% Ideally you will provide here the full path to the saving directory.
% Something like:
% saving_directory = '/Users/Anush/GIT_Folders/phenomapping/tmpresults/';
% Alternatively, since you work from the parent phenomapping folder you can
% define:
% saving_directory = 'tmpresults/';


% See in the comment below an example of the paths to phenomapping, matTFA
% and TEX-FBA repositories and CPLEX
% phenomapping_directory = '/Users/Anush/GIT_Folders/phenomapping';
% mattfa_directory = '/Users/Anush/GIT_Folders/matTFA';
% texfba_directory = '/Users/Anush/GIT_Folders/texfba';
% cplex_directory = '/Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271';

% Add path to phenomapping repository:
% addpath(genpath(phenomapping_directory));

% Call the function "initPhenoMappingPaths.m" to add all other required 
% paths
% [phenomapping_directory, thermo_data_directory] = ...
% initPhenoMappingPaths(saving_directory);

% You will also need to paste the path to cplex when the following message
% appears:
% "Please provide your cplex path and press enter"
% You will paste the corresponding path as follows (without " "):
% ... /Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271

% Work from parent folder of phenomapping as the starting path. If you are
% not there run the following command:
% cd(phenomapping_directory)

% Tag true to apply phenomapping with TFA (thermodynamic constraints will
% be taken into account)
% tagThermo = 1;

% Upload model - here iPbe for liver-stage analysis provided in the
% subfolder models as test case
% modeldescription = 'iPbeLiver';
% modelPath = '/Users/Anush/GIT_Folders/phenomapping/models/tipbe4liver.mat';
% load(modelPath)
% model = tipbe;

% Prepare model for phenomapping
% [model, checkList, tagReady] = initTestPhenoMappingModel(model,...
%     thermo_data_directory,tagThermo);

% Save prepared model
% save(strcat(saving_directory,'inputPhenoMapping.mat'), 'model');

%% PhenoMapping analysis for condition specific information in the GEM
% as applied for the analysis of iPbe and construction of iPbe-liver and 
% iPbe-blood

% Work from parent folder of phenomapping as the starting path. If you are
% not there run the following command:
% cd(phenomapping_directory)

% Step 0: Get essentiality for the non-stage specific or generic model. 
% This should be a model without context-specific data integrated. 
% iPbe is an example of such a generic model. 
% iPbe is saved in the files tipbe4blood.m and tipbe4liver.m. Both versions
% are stoichiometrically identical. They differ in the thermodynamic data
% integrated, i.e., they integrate pH and membrane potentials for the blood- 
% and liver-stages of the malaria infection. This allows that the 
% consequent phenomapping analyses are thermodynamically consistent with
% the corresponding life stage.

% if you are running the core substrates directly after settings
% you don't need to uncomment this. Else uncomment:

% modeldescription = 'iPbeLiver';
% filename = strcat(modeldescription,'_PhenoMappingSettings');
% saving_directory = strrep(which(strcat(filename, '.mat')), ...
%     strcat(filename, '.mat'), '');
% clearvars -except filename saving_directory modeldescription
% load(strcat(saving_directory,filename,'.mat'));

run('tests/test_core_essentiality.m')


% Step 1: Map essentiality to substrate availability following the 
% PhenoMapping analysis for the medium:

% if you are running the core substrates directly after core essentiality
% you don't need to uncomment this. Else uncomment:

% modeldescription = 'iPbeLiver';
% filename = strcat(modeldescription,'_PhenoMappingEssentiality');
% saving_directory = strrep(which(strcat(filename, '.mat')), ...
%     strcat(filename, '.mat'), '');
% clearvars -except filename saving_directory
% load(strcat(saving_directory,filename,'.mat'));

run('tests/test_core_substrates.m')


% Step 2: Map essentiality to metabolite concentrations following the 
% PhenoMapping analysis for metabolomics (if available):

modeldescription = 'iPbeLiver';
filename = strcat(modeldescription,'_PhenoMappingEssentiality');
saving_directory = strrep(which(strcat(filename, '.mat')), ...
    strcat(filename, '.mat'), '');
clearvars -except filename saving_directory
load(strcat(saving_directory,filename,'.mat'));

metabolomics_description = {'metabolomics blood pfa'};
metabolomics_directory = which('allmetab_pfa_blood.mat');

run('tests/test_core_metabolomics.m')


% Step 3: Map essentiality to gene expression levels following the 
% PhenoMapping analysis for transcriptomics (if available) considering or 
% not regulation of gene expression between isoenzymes:
% We perform Thermodynamic-based flux Variability Analysis (TVA) to get the
% reaction flux ranges that will be an input to the TEX-FBA
modeldescription = 'iPbeLiver';
filename = strcat(modeldescription,'_PhenoMappingEssentiality');
saving_directory = strrep(which(strcat(filename, '.mat')), ...
    strcat(filename, '.mat'), '');
clearvars -except filename saving_directory
load(strcat(saving_directory,filename,'.mat'));

transcriptomics_description = {'transcriptomics liver pbe'};
transcriptomics_directory = which(...
    'levelGenes_pbe_liver_HepG2_mean_48h.mat');

run('tests/test_core_mattfa_minmax.m')
run('tests/test_core_transcriptomics.m')

%% Extract data stored in tmpresults
% The files stores in mat structures in tmpresults can be read an the 
% information extracted into text or csv files for further manual analysis
modeldescription = 'iPbeLiver';
filename = strcat(modeldescription,'_PhenoMappingSettings');
saving_directory = strrep(which(strcat(filename, '.mat')), ...
    strcat(filename, '.mat'), '');
clearvars -except filename saving_directory modeldescription
load(strcat(saving_directory,filename,'.mat'));

run('tests/test_io.m')
