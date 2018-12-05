function [DPs, model, objectives] = findDPMinMets(model, NumAlt, indUSE, time, tagMin, filename)
% Get alternative solutions for MILP (maximization)
%
% USAGE:
%
%       [DPs, model, objectives] = findDPMinMets(model, NumAlt, indUSE, time, tagMin, filename)
%
% INPUTS:
%    model:           Model with TFA structure and MILP formulation
%
% OPTIONAL INPUTS:
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    indUSE:          Indexes of integers in the MILP (default = 'BFUSE')
%    time:            Time in sec after which simulation will be stopped
%                     (default = empty / no time provided)
%    tagMin:          True to generate alternative solution of minimal size
%                     only. False to generate up to the number of
%                     alternatives provided in NumAlt (default = true)
%    filename:        Name used to save DPs (default = 'imm')
%                     
%
% OUTPUTS:
%    DPs:             Directionality profile matrix with alternatives in
%                     each column
%    model:           Model with integer cuts integrated to avoid
%                     repetition of same solution or supersolution (this 
%                     is a superset of an obtained solution)
%    objectives:      Exchange/Drain reactions
%
% .. Author:
% Anush Chiappino 2017
% 

if (nargin < 2)
    NumAlt = 1;
end
if (nargin < 3) || isempty(indUSE)
    indUSE = model.indUSE;
end
if (nargin < 4)
    time = [];
end
if (nargin < 5)
    tagMin = 1;
end
if (nargin < 6)
    filename = 'imm';
end

NumSols = 0;
sol = optimizeThermoModel(model,0,'cplex',time);
DPs = [];
objectives = [];

if ~(isempty(sol.x)) && tagMin
    [NumCons,NumVars] = size(model.A);
    NewCons = zeros(NumVars,1);
    NewCons(indUSE) = 1;
    model.A(NumCons+1,:) = NewCons;
    model.rhs(NumCons+1) = sol.val-0.5;
    model.constraintNames{NumCons+1} = ['CUT_0'];
    model.constraintType{NumCons+1} = '>';
end

while ((NumSols < NumAlt) && ~(isempty(sol.x)) && ~(sol.val==0))
    [NumCons,NumVars] = size(model.A);
    
    if ~(isempty(sol.x))
        NumSols = NumSols + 1;
        objectives(NumSols,1) = sol.val;
        DPs(:,NumSols) = sol.x;
        
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
        sol = optimizeThermoModel(model,0,'cplex',time);
        
        if isempty(sol.x) || (sol.val==0)
            break;
        end
        fprintf('Number of DPs:\t%d\n',NumSols);
        
        if rem(NumSols,5) == 0
            save(strcat('tmpresults/',filename,'_DPs.mat'), 'DPs');
        end
    end
end
if (NumSols == NumAlt)
    sol = optimizeThermoModel(model,0,'cplex',time);
    if ~isempty(sol.x) || (sol.val>0)
        warning('You need to generate more alternative solutions. Check tutorial_issues.m');
    end
end
end




