function model = recoverModel4DPMax(modelmilp, DPs)
% Integrate integer cuts into model from DPs matrix (needed if matlab
% crashes) for a maximization problem
%
% USAGE:
%
%    model = recoverModel4DP(modelmilp, DPs)
%
% INPUT:
%    model4DP:    TFA model with MILP structure for a maximization problem
%    DPs:         Directionality profile matrix with alternatives in
%                 each column from the problem of study 
%
% OUTPUTS:
%    model:       TFA model with MILP structure for a maximization problem 
%                 and integet cuts integrated
%
% .. Author:
% Meric Ataman & Anush Chiappino-Pepe 2017
% 

model = modelmilp;
intTag = {'BFUSE'};
indUSE = getAllVar(model,intTag);

for NumSols = 1:size(DPs,2)
    sol.x = DPs(:,NumSols);
    [NumCons,NumVars] = size(model.A);
    
    % we find all the use vectors and formulate them into a new integer cut
    % constraint
    USEvecSol = ones(NumVars,1);
    USEvecSol(indUSE) = sol.x(indUSE);
    actUSEvec = find(USEvecSol<0.1);
    
    cardinality = length(actUSEvec);
    NewCons = zeros(NumVars,1);
    NewCons(actUSEvec) = 1;
    
    model.A(NumCons+1,:) = NewCons;
    model.rhs(NumCons+1) = 0.5;
    model.constraintNames{NumCons+1} = ['CUT_' num2str(NumSols)];
    model.constraintType{NumCons+1} = '>';
    
end
end
            