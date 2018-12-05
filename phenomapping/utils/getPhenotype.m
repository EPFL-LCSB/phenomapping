function phenotype = getPhenotype(listGenes,lifeStage,PlasmoSpecie,pathToData)
% Get phenotypes for a list of genes and path provided (default: gets
% Plasmodium berghei blood and liver stage phenotypes)
%
% USAGE:
%
%    phenotype = getPhenotype(listGenes,lifeStage,PlasmoSpecie,pathToData)
%
% INPUT:
%    listGenes:       cell array of geneIDs
%
% OPTIONAL INPUTS:
%    lifeStage:       Either 'blood' or 'liver' (Default = 'blood')
%    PlasmoSpecie:    Either 'pbe' or 'pfa' (Default = 'pbe')
%    pathToData:      Cell array with path(s) to phenotype(s) for a life 
%                     stage of the organism of study (default = empty / get
%                     data for P. berghei blood and liver stages saved in 
%                     phenomapping/data/pbe). NOTE: new phenotypes should
%                     be saved in a cell called "phenotypes".
%
% OUTPUTS:
%    phenotype:       Phenotypes extracted
%
% .. Author:
%       - Anush Chiappino-Pepe 2017

if (nargin < 2)
    lifeStage = 'blood';
end
if (nargin < 3)
    PlasmoSpecie = 'pbe';
end
if (nargin < 4)
    pathToData = [];
end

listGenes = strrep(listGenes,'putative_','');
listGenes = strtok(listGenes,'.');

if isempty(pathToData)
    if strcmp(PlasmoSpecie,'pbe')
        if strcmp(lifeStage,'blood')
            mm = load('/phenomapping/data/pbe/pbe_phenotypes_blood_Nov16.mat');
            spPheno = mm.phenotypes;
            clear mm
        elseif strcmp(lifeStage,'liver')
            mm = load('/phenomapping/data/pbe/pbe_phenotypes_liver_Jul18.mat');
            spPheno = mm.phenotypes;
            clear mm
        end
    end
else
    mm = load(pathToData);
    spPheno = mm.phenotypes;
    clear mm
end

phenotype = cell(length(listGenes), size(spPheno,2));
[y,r]=ismember(listGenes,spPheno(:,1));
phenotype(y,:) = spPheno(r(y),:);
phenotype(~y,:) = {'no info'};
end
