function data = exportPhenoMappingInfoTable(essGenes,conditionToGene,conditionDesc)
% Generates table of output of PhenoMapping (linkEssGeneIMM2Subs.m,
% getBotNeckMets)
%
%
% USAGE:
%
%    data = exportPhenoMappingInfo(essGenes,conditionToGene,conditionDesc)
%
% INPUT:
%    essGenes:        List of genes that are essential in a condition and
%                     not in the reference model
%    conditionToGene: Conditions underlying essential genes based on
%                     phenomapping
%
% OPTIONAL INPUTS:
%    conditionDesc:   Cell array with description of condition (default
%                      = 'condition in PhenoMapping')
%
%
% OUTPUTS:
%    table:           Summary
%
% .. Author:
%       - Anush Chiappino-Pepe 2019
%

if (nargin < 3) || isempty(conditionDesc)
    conditionDesc = {'condition in PhenoMapping'};
end

[conditionToGene] = collectInfoCell(conditionToGene);

% identify max number of alternatives found in conditionToGene
NumAlt = 1;
for i = 1:size(conditionToGene,1)
    NumAlt = max([NumAlt;size(conditionToGene{i},2)]);
end

% create heading for alternatives
heading = cell(1,NumAlt);
for i = 1:NumAlt
    heading(i) = strcat('Alt',{' '},num2str(i));
end

data = [];
for i = 1:size(conditionToGene,1)
    if ~isempty(conditionToGene{i})
        for j = 1:size(conditionToGene{i},1)
            for z = 1:NumAlt
                if z>(size(conditionToGene{i},2)) || ...
                        isempty(conditionToGene{i}{j,z})
                    conditionToGene{i}{j,z} = '';
                else
                    conditionToGene{i}{j,z} = strrep(conditionToGene{i}...
                        {j,z},strcat({' '},'<=>',{' '}),'');
                    conditionToGene{i}{j,z} = strrep(conditionToGene{i}...
                        {j,z},strcat({' '},'<=>'),'');
                    conditionToGene{i}{j,z} = strrep(conditionToGene{i}...
                        {j,z},',','');
                end
                if iscell(conditionToGene{i}{j,z})
                    conditionToGene{i}{j,z} = conditionToGene{i}{j,z}{1};
                end
            end
        end
        data = [data; [strcat(conditionDesc,{' '},num2str(i)), ...
            heading]; essGenes{i},conditionToGene{i}];
    end
end
