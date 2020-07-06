function exportPhenoMappingInfo(essGenes,conditionToGene,conditionDesc,pathToPhenoData,phenDesc,filename,extension,tagPhenotype)
% Exports output of PhenoMapping (linkEssGeneIMM2Subs.m, getBotNeckMets)
% to text or cvs files and links genes essential at the conditions with 
% available phenotypes
% 
%
% USAGE:
%
%    exportPhenoMappingInfo(essGenes,conditionToGene,conditionDesc,pathToPhenoData,phenDesc,filename)
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
%    pathToPhenoData: Cell array with path(s) to phenotype(s) data for the 
%                     organism of study (default = empty / provide a path). 
%                     NOTE: the phenotypes should be saved in a cell called
%                     "phenotypes": 1st and 2nd column geneIDs and 
%                     phenotypes, respectively.
%    phenDesc:        Cell array with description of phenotype(s) (default
%                      = 'phenotype')
%    filename:        Name used to save file (default = 
%                     'PhenoMappingCondition2Genes'). NOTE: no need to
%                     provide the extension of the file here. It will be
%                     csv by default. But the path where the file will be
%                     saved should be provided else it will be saved in the
%                     folder where you currently are.
%    extension:       Text (txt) or csv (csv) file (default = 'txt')
%    tagPhenotype:    True to map phenotype (default)
%                     
%
% OUTPUTS:
%    cvs/text:        Saved in phenomapping/tmpresults. You can open the  
%                     file in excel and separate Text to columns with Tab
%                     as a delimiter.
%
% .. Author:
%       - Anush Chiappino-Pepe 2018
%

if (nargin < 3) || isempty(conditionDesc)
    conditionDesc = {'condition in PhenoMapping'};
end
if (nargin < 4) || all(strcmp(pathToPhenoData,''))
    pathToPhenoData = [];
end
if (nargin < 5) || all(isempty(phenDesc)) || all(strcmp(phenDesc,''))
    phenDesc = pathToPhenoData;
end
if (nargin < 6) || isempty(filename)
    filename = 'PhenoMappingCondition2Genes';
end
if (nargin < 7) || isempty(extension)
    extension = 'txt';
end
if (nargin < 8) || isempty(tagPhenotype)
    tagPhenotype = 1;
end

% converting characters to cells to avoid issues
if ischar(pathToPhenoData)
    pathToPhenoData = {pathToPhenoData};
end

if (size(phenDesc,1)>size(phenDesc,2))
    phenDesc = phenDesc';
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
        if tagPhenotype
            if ~isempty(pathToPhenoData)
                phen = [];
                for j = 1:length(pathToPhenoData)
                    tmp = getPhenotype(essGenes{i},'','',pathToPhenoData{j});
                    phen = [phen, tmp(:,2)];
                end
            else
                phenDesc = {'phenotype'};
                tmp = getPhenotype(essGenes{i},'','',[]);
                phen = tmp(:,2);
            end
            data = [data; [strcat(conditionDesc,{' '},num2str(i)),phenDesc, ...
                heading]; essGenes{i},phen,conditionToGene{i}];
        else
            data = [data; [strcat(conditionDesc,{' '},num2str(i)), ...
                heading]; essGenes{i},conditionToGene{i}];
        end
    end
end


desc = '%s\t';
for i = 2:size(data,2)
    desc = strcat(desc,'%s\t');
end

writeData(strcat(filename,'.',extension),data,desc);