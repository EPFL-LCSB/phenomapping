% You can find here SOPs on how to deal with some issues that might arise during the IMM analysis

%% if matlab crashes while generating the alternative IMMs:
% you dont need to regenerate all alternatives again 
% you can recover the cut constraints from the DPs (saved in tmpresults) 
% and the modelmilp (output from analysisIMM)
if exist(strcat('tmpresults/',filename,'_DPs.mat'),'file') == 2
    DPr = load(strcat('tmpresults/',filename,'_DPs.mat'));
    model4DP = recoverModel4DPIMM(modelmilp, DPr.DPs);
else
    fprintf('no DPs-matrix for IMM was found\n');
end

%% if you need to generate more alternative IMMs:
% you need to call again the function findDPMinMets providing a high number
% of alternatives (suggested NumAlt = 5000)
% you should use here the model4DP from the previous generation of 
% findDPMinMets - this will allow to integrate new cut constraints on the
% top of the previous ones
% to consider: for models with multiple nutrient possibilities 
% (like parasites) and/or that are very relaxed, the number of alternative
% IMMs can go up to 10000s
[DPs, model4DP] = findDPMinMets(model4DP, NumAlt, modelmilp.indUSE, ...
       time, filename);
   
% Note: the input and output notation of this script is consistent with
% the notation of the Example_script_imm.m

