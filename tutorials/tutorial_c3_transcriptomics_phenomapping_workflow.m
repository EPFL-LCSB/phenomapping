%% Bottleneck genes alternatives

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/analysis_blood/bottleneckgenes/ipfav2red_bloodcsallmets_ge_input_botgenes.mat')
% this file contains the expModelm, solution (from ipfav2red_bloodcsallmets_ge_recon), and addEssGenes
clear
clc

% mm = load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbev2/analyses_liver/gene_exp_liver/ipbe2red_liver_ge_output_ess.mat');
mm = load('/Users/Anush/Dropbox/Aarti&Anush/analysis_itgo_Sep18/analysis/gene_exp/itgo_brady_ge_ess.mat');
expModelm = mm.expModelm;
solution = mm.solution;
tagReg = 0;

% mmorig = load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbev2/analyses_liver/essentiality/ess_tipbe2_liver.mat');
mmorig = load('/Users/Anush/Dropbox/Aarti&Anush/analysis_itgo_Sep18/analysis/essentiality/itgo_ess.mat');
essorig = mmorig.essential_genes_tfa;
clear mmorig

% should be the same length as the consistency score
grRate = 0.01;
indNF = getAllVar(expModelm,{'NF'});
numAlt = 20;
indUPDOWN = [expModelm.ExpInfo.indUP; expModelm.ExpInfo.indDOWN];
essThr = 0.1;

if tagReg
    col = 1;
else
    col = 3;
end

allToRelax = cell(1,size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2));
addEssGenesAs = cell(1,size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2));

for selectAlti = 1:size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2) % select alternative to integrate
    newess = mm.essGenes_expModelAlt{selectAlti,col};
    addEssGenes = newess(not(ismember(newess,essorig))); %additional genes identified with gene expression data integrated vs previous version (TFA tipbe2 liver)
    addEssGenesAs{1,selectAlti} = addEssGenes;    
    indObj = indUPDOWN(solution.sol_matrix(indUPDOWN,selectAlti)>0.9);
    for i = 1:length(addEssGenes)
        modeli = thermoDeleteModelGenes(expModelm,addEssGenes(i));
%         [~, ~, ~, ~, modeli] =
%         thermoSingleGeneDeletionGeneExpNewGPR(expModelm, grRate, essThr,
%         0, addEssGenes(i), 0, indObj,numAlt,indNF); % this analysis
%         doesnt work for genes identified exclusively with the NOREG
%         analysis
        [genesToKeep,genesToRelax] = getBotGenes(modeli,indObj,grRate,numAlt,indNF);
        allToRelax{i,selectAlti} = genesToRelax;
    end
end

clear modeli
clear genesToKeep
clear genesToRelax
% [genesToKeep,genesToRelax] = getBotGenes(expModel,indObj,grRate,numAlt,indNF);

allToRelaxLump = cell(size(allToRelax,1),size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2));
for selectAlti = 1:size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2) % select alternative to integrate
    for i = 1:size(allToRelax,1)
        for j = 1:length(allToRelax{i,selectAlti})
            if isempty(allToRelax{i,selectAlti}{j})
                allToRelaxLump{i,selectAlti}{j} = '';
            else
                allToRelaxLump{i,selectAlti}{j} = concatenateList(allToRelax{i,selectAlti}{j}, '|');
            end
        end
        if length(allToRelaxLump{i,selectAlti})>1.5
            allToRelaxLump{i,selectAlti} = concatenateList(allToRelaxLump{i,selectAlti}, ':NEXT-ALT:');
        end
    end
end

uniqueGenes = [];
for i = 1:size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2)
    uniqueGenes = unique([uniqueGenes; addEssGenesAs{1,i}]);
end

sumAllToRelaxLump = cell(size(uniqueGenes,1),size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2));
for i = 1:size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2)
    [y,r] = ismember(uniqueGenes,addEssGenesAs{1,i});
    sumAllToRelaxLump(y,i) = allToRelaxLump(r(y),i);
end

sumAllToRelaxLumpSame = zeros(size(uniqueGenes,1),1);
for j = 1:size(uniqueGenes,1)
    y = zeros(size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2)-1,1);
    for i = 2:size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2)
        y(i) = strcmp(sumAllToRelaxLump{j,i},sumAllToRelaxLump{j,1});
    end
    if sum(y)==size(solution.z_matrix(solution.store_obj>max(solution.store_obj)-0.5),2)-1
        sumAllToRelaxLumpSame(j) = 1;
    end
end
