function [phenomapping_directory, thermo_data_directory] = ...
    initPhenoMappingPaths(saving_directory, mattfa_directory, ...
    texfba_directory)
% Initialization of paths for phenomapping
%
% USAGE:
%
%    [phenomapping_directory, thermo_data_directory] = initPhenoMappingPaths(saving_directory)
%
% OPTIONAL INPUTS:
%    saving_directory:      Directory where we want to save intermediate  
%                           and final results (default = tmpresults folder  
%                           in parent phenomapping folder)
%    mattfa_directory:      Directory to the matTFA repository (default = 
%                           empty / provide manually)
%    texfba_directory:      Directory to the texfba repository (default = 
%                           empty / provide manually)
%
% OUTPUTS:
%    texfba_directory:      Directory to the parent phenomapping folder
%    thermo_data_directory: Directory to the thermodynamic data within
%                           matTFA
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

if (nargin < 1)
    saving_directory = strrep(mfilename('fullpath'),...
    'phenomapping/initPhenoMappingPaths','tmpresults/');
end
if (nargin < 2)
    mattfa_directory = [];
end
if (nargin < 3)
    texfba_directory = [];
end

phenomapping_directory = strrep(mfilename('fullpath'),...
    'phenomapping/initPhenoMappingPaths','');
addpath(genpath(phenomapping_directory));

if isempty(mattfa_directory)
    mattfa_directory = uigetdir('matTFA','select the matTFA directory');
end
addpath(genpath(mattfa_directory));
thermo_data_directory = strcat(mattfa_directory,...
    '/thermoDatabases/thermo_data.mat');

if isempty(texfba_directory)
    texfba_directory = uigetdir('texfba','select the texfba directory');
end
addpath(genpath(texfba_directory));

cd(phenomapping_directory)

% PhenoMapping was developed to work with the solver CPLEX. We hence check 
% that you have CPLEX installed. In future releases, this repository
% will work with other solvers like gurobi.
cplex_directory = what('cplex');
if isempty(cplex_directory)
    cplex_directory = [];
else
    cplex_directory = cplex_directory.path;
end
addCplexPath(cplex_directory);

% create temporary results folder to save intermediate and final results: 
% note that PhenoMapping will not upload the intermediate results 
% - but you might need to use them if matlab crashes
if ~isdir(saving_directory(1:end-1))
    mkdir(saving_directory(1:end-1))
end