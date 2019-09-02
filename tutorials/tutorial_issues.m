% You can find here SOPs on how to deal with some issues that might arise 
% during the PhenoMapping analyses, which involves several Mixed Integer 
% Linear Programming (MILP) formulations

%% if matlab crashes while generating the alternative MILP solutions:
% - You dont need to regenerate all alternative solutions again.
% 
% - You can recover solutions generated up to that point and the 
% corresponding cut constraints, which basically define that the already 
% identified solutions will not be found again.
%
% - This is done by first uploading the directionality profiles or DPs  
% saved in the subfolder "tmpresults", and the model with the milp  
% formulation, e.g. modelmilp (output from analysisIMM). See below an
% example

if exist(strcat('tmpresults/',filename,'_DPs.mat'),'file') == 2
    DPr = load(strcat('tmpresults/',filename,'_DPs.mat'));
    model4DP = recoverModel4DPMax(modelmilp, DPr.DPsimm);
else
    fprintf('no DPs-matrix for IMM was found\n');
end

% - Call then the function "recoverModel4DPMax.m" as suggested
% below

[DPsimm, model4DP2] = findDPMax(model4DP, NumAlt, modelmilp.indUSE, ...
       time, 1, filename);
   
% Note: the input and output notation provided here is consistent with
% the notation of the test_core_substrates.m
