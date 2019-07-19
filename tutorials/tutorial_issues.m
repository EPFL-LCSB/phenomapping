% You can find here SOPs on how to deal with some issues that might arise 
% during the Mixed Integer Linear Programming (MILP) analysis

%% if matlab crashes while generating the alternative MILP solutions:
% you dont need to regenerate all alternative solutions again 
% you can recover the cut constraints from the DPs (saved in tmpresults) 
% and the model with the milp formulation e.g. modelmilp (output from 
% analysisIMM)
% For that purpose use the function %%recoverModel4DPMax%% as suggested
% below
if exist(strcat('tmpresults/',filename,'_DPs.mat'),'file') == 2
    DPr = load(strcat('tmpresults/',filename,'_DPs.mat'));
    model4DP = recoverModel4DPMax(modelmilp, DPr.DPsimm);
else
    fprintf('no DPs-matrix for IMM was found\n');
end

%% if you need to generate more alternative MILP solutions:
% you need to call again the function findDPMax or findDPMin providing a 
% high number of alternatives (suggested NumAlt = 5000)
% you should use here the model from the previous generation of 
% findDPMax or findDPMin - this will allow to integrate new cut constraints 
% on the top of the previous ones
% to consider: the number of MILP alternatives can be very high (orders of
% magnitude), especially if the model is not very constrained and the
% number of integer variables is high
[DPsimm, model4DP2] = findDPMax(model4DP, NumAlt, modelmilp.indUSE, ...
       time, 1, filename);
   
% Note: the input and output notation of this script is consistent with
% the notation of the test_core_substrates.m

