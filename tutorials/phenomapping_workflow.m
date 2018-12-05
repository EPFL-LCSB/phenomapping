%% Setting up the solver and preparing the model for PhenoMapping
% Step 1: PhenoMapping was developed to work with the solver CPLEX. As a 
% first step we hence check that you have CPLEX installed. In future
% releases, this repository will work with other solvers like gurobi.
cplexPath = '/Users/Anush/Applications/IBM/ILOG/CPLEX_Studio1271';
solverOK = addCplexPath(cplexPath);

% Step 2: check the structure of your model and prepare it for a
% PhenoMapping analysis
modelPath = '';
load(modelPath)
modeli = modelname;
[model, checkList, tagReady] = initTestPhenoMappingModel(modeli);

%% PhenoMapping analysis for organism specific information in the GEM
% Step 1: Map essentiality to enzymatic irreversibilities defined as adhoc 
% in the model (if applicable): 
% b1_adhocirrev_phenomapping_workflow.m


% Step 2: Map essentiality to limitation in intracellular transportability
% (if eukaryotic): 
% b2_inttransport_phenomapping_workflow.m


% Step 3: Map essentiality to localization of enzymes (if eukaryotic): 
% b3_enyloc_phenomapping_workflow.m


% Step 4: Map essentiality to biochemistry annotated to genome: 
% b4_biochem_phenomapping_workflow.m



%% PhenoMapping analysis for condition specific information in the GEM
% Step 1: Map essentiality to substrate availability following the 
% PhenoMapping analysis for the medium: 
c1_medium_phenomapping_workflow.m


% Step 2: Map essentiality to metabolite concentrations following the 
% PhenoMapping analysis for metabolomics (if available):
% c2_metabolomics_phenomapping_workflow.m


% Step 3: Map essentiality to gene expression levels following the 
% PhenoMapping analysis for transcriptomics (if available) considering or 
% not regulation of gene expression between isoenzymes: 
% c3_transcriptomics_phenomapping_workflow.m

