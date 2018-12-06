function phenotype = getPhenotype(listGenes,lifeStage,PlasmoSpecie,pathToData)
% Get phenotypes for a list of genes and path provided (default: provide 
% path to phenotypes). This allows to save a path and upload files based 
% on the id of an organism and the context of study. It was implemented for
% Plasmodium berghei (pbe) blood (blood) and liver (liver) stages.
%
% USAGE:
%
%    phenotype = getPhenotype(listGenes,lifeStage,PlasmoSpecie,pathToData)
%
% INPUT:
%    listGenes:       cell array of geneIDs
%
% OPTIONAL INPUTS:
%    lifeStage:       Either 'blood' or 'liver' (Default = '')
%    PlasmoSpecie:    Either 'pbe' or 'pfa' (Default = '')
%    pathToData:      Cell array with path(s) to phenotype(s) for a life
%                     stage of the organism of study (default = empty / get
%                     data for P. berghei blood and liver stages saved in
%                     phenomapping/data/pbe). NOTE: new phenotypes should
%                     be saved in a cell called "phenotypes".
%
% OUTPUTS:
%    phenotype:       Phenotypes extracted for the list of genes provided
%
% .. Author:
%       - Anush Chiappino-Pepe 2017

if (nargin < 2)
    lifeStage = '';
end
if (nargin < 3)
    PlasmoSpecie = '';
end
if (nargin < 4)
    pathToData = [];
end

if ischar(listGenes)
    listGenes = {listGenes};
end

curdir = cd;

listGenes = strrep(listGenes,'putative_','');
listGenes = strtok(listGenes,'.');

if strcmp(PlasmoSpecie,'pbe')
    if strcmp(lifeStage,'blood')
        pathToData = strcat(curdir,'/data/pbe/pbe_phenotypes_blood_Nov16.mat');
        if ~exist(pathToData,'file') == 2
            pathToData = [];
        end
        while isempty(pathToData)
            pathToData = input('Please provide the path to the blood phenotypes of P. berghei and press enter\n... ','s');
        end
    elseif strcmp(lifeStage,'liver')
        pathToData = strcat(curdir,'/data/pbe/pbe_phenotypes_liver_Jul18.mat');
        if ~exist(pathToData,'file') == 2
            pathToData = [];
        end
        while isempty(pathToData)
            pathToData = input('Please provide the path to the liver phenotypes of P. berghei and press enter\n... ','s');
        end
    end
else
    while isempty(pathToData)
        pathToData = input('Please provide the path to the phenotypes of interest and press enter\n... ','s');
    end
end

mm = load(pathToData);
if ~isfield(mm,'phenotypes')
    error('the phenotypes should be saved in a two-column cell called "phenotypes": the first column are the geneIDs and the second column the phenotypes')
end
spPheno = mm.phenotypes;
clear mm

phenotype = cell(length(listGenes), size(spPheno,2));
[y,r] = ismember(listGenes,spPheno(:,1));
phenotype(y,:) = spPheno(r(y),:);
phenotype(~y,:) = {'no info'};
end
