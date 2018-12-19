function [modeldescription, saving_directory] = retrievePaths...
    (modeldescription, saving_directory)
% Short function to retrieve saving path and model description. Just to
% make sure the user provides this information before it is required.
%
% USAGE:
%
%    [modeldescription, saving_directory] = retrievePaths(modeldescription, saving_directory)
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

if ~exist('modeldescription','var') || nargin < 1
    modeldescription = [];
end
if ~exist('saving_directory','var') || nargin < 2
    saving_directory = [];
end
while isempty(modeldescription)
    modeldescription = input('Please upload PhenoMappingSettings or provide modeldescription and press enter\n... ','s');
end
while isempty(saving_directory)
    saving_directory = input('Please upload PhenoMappingSettings or provide saving_directory and press enter\n... ','s');
end

