function [DPs, model, objectives] = findDPMin(model, NumAlt, indUSE, time, tagMax, filename)
% Get alternative solutions for MILP (minimization)
%
% USAGE:
%
%       [DPs, model, objectives] = findDPMin(model, NumAlt, indUSE, time, tagMax, filename)
%
% INPUTS:
%    model:           Model with TFA structure and MILP formulation
%
% OPTIONAL INPUTS:
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    indUSE:          Indexes of integers in the MILP (default = 
%                     model.indUSE)
%    time:            Time in sec after which simulation will be stopped
%                     (default = empty / no time provided)
%    tagMax:          True to generate alternative solution of minimal size
%                     only. False to generate up to the number of
%                     alternatives provided in NumAlt (default = false)
%    filename:        Name used to save DPs (default = 'PhenoMappingDPMin')
%                     
%
% OUTPUTS:
%    DPs:             Directionality profile matrix with alternatives in
%                     each column
%    model:           Model with integer cuts integrated to avoid
%                     repetition of same solution or subset of the solution
%    objectives:      Optimal value, which is the sum of active intigers,
%                     after the optimization
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
    tagMax = 0;
end
if (nargin < 6)
    filename = 'PhenoMappingDPMin';
end

[~,NumVars] = size(model.A);
model.f = zeros(NumVars,1);
model.f(indUSE) = 1;

model.objtype = 1; % minimization
NumSols = 0;
sol = optimizeThermoModel(model,0,'cplex',time);
DPs = [];
objectives = [];

if ~(isempty(sol.x)) && tagMax
    [NumCons,NumVars] = size(model.A);
    NewCons = zeros(NumVars,1);
    NewCons(indUSE) = 1;
    model.A(NumCons+1,:) = NewCons;
    model.rhs(NumCons+1) = sol.val+0.5;
    model.constraintNames{NumCons+1} = ['CUT_0'];
    model.constraintType{NumCons+1} = '<';
end

while ((NumSols < NumAlt) && ~(isempty(sol.x)))
    
    [NumCons,NumVars] = size(model.A);
    
    if ~(isempty(sol.x))
        NumSols = NumSols + 1;
        objectives(NumSols,1) = sol.val;
        DPs(:,NumSols) = sol.x;
        USEvecSol = zeros(NumVars,1);
        USEvecSol(indUSE) = sol.x(indUSE);
        actUSEvec = find(USEvecSol > 0.98);
        
        cardinality = length(actUSEvec);
        NewCons = zeros(NumVars,1);
        NewCons(actUSEvec) = 1;
        
        model.A(NumCons+1,:) = NewCons;
        model.rhs(NumCons+1) = cardinality - 0.5;
        model.constraintNames{NumCons+1} = ['CUT_' num2str(NumSols)];
        model.constraintType{NumCons+1} = '<';
        
        % solve the problem again
        sol = optimizeThermoModel(model,0,'cplex',time);
        
        fprintf('Number of DPs:\t%d\n',NumSols);
        
        if rem(NumSols,10) == 0
            save(strcat(filename,'_DPs.mat'), 'DPs');
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




