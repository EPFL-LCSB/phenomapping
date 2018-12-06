clear
clc

%% load model
% iPbe
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbev2/model_liver/tipbe2_liver.mat')
% model = tipbe_liver;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbev2/model_liver/tipbe2_liver_reducedModel.mat')
% model = treducedModel;

% description = 'iPbe2 liver cs';
% mm = load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbev2/model_liver/tipbe2_liver_cs_reducedModel.mat');
% model = mm.treducedModel;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbe/model_blood/tipbe_blood.mat')
% model = tipbe;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbe/model_blood/treducedModel_blood.mat')
% model = treducedModel;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPbe/model_blood/treducedModel_bloodcs_allmets.mat')
% model = treducedModel_bcsallmets;


% iPfav2
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/model_liver/tipfav2_liver_corr.mat')
% model = tipfav2;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/model_liver/tipfav2_liver_reducedModel_corr.mat')
% model = treducedModel;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/model_blood/tipfav2_blood_corr.mat')
% model = tipfav2;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/model_blood/tipfav2_blood_reducedModel_corr.mat')
% model = treducedModel;

% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/model_blood/treducedModel_bloodcs_allmets.mat')
% model = treducedModel_bcsallmets;

description = 'iTgo tachy';
% description = 'iTgo brady';
mm = load('/Users/Anush/Dropbox/Aarti&Anush/analysis_itgo_Sep18/itgo_reducedModel.mat');
model = mm.treducedModel;

%% load RNAseq data

% blood
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/14_pbe/data/Pbe_blood_RNAseq.mat')
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/data/2010_Otto_SI4_pfa_data.mat')

% liver A
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/14_pbe/data/Counts_OnPbANKA_2a_Sunil_Feb17.mat')
% liver B
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/14_pbe/data/Counts_OnPbANKA_2b_Sunil_Feb17.mat')
% genesrnainfo=sunil; % comment for blood
% genesrnainfovalue=sunilval;

% liver pfa

% hypnozoites
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/dormancy/data/bhatia_data_topfa.mat')
% genesrnainfo = genesfal;
% genesrnainfo=strtok(genesrnainfo,'.');
% genesrnainfo=strrep(genesrnainfo,{' '},'');
% genesrnainfovalue = numfal;
% level = 1;


load('/Users/Anush/Dropbox/Aarti&Anush/analysis_itgo_Sep18/data/RNAseq_data_toxo.mat')
genesrnainfo = genes;
genesrnainfo=strtok(genesrnainfo,'.');
genesrnainfo=strrep(genesrnainfo,{' '},'');

% tachy Toxo Aarti
genesrnainfovalue = values(:,1);
level = 1;

% brady Toxo Aarti
% genesrnainfovalue = values(:,2);
% level = 1;


%% decide time point
% level=3; %3:4=trophozoite Pbe in blood, 4 trophozoite pfa blood, 3=48h in liver

%% preprocessing of genes in model

genes=strrep(model.genes,'putative_','');
genes=strtok(genes,'.');
genes=strrep(genes,{' '},'');

% blood and liver
[yes,row]=ismember(genes,genesrnainfo);
sum(yes)
levelGenes=zeros(length(model.genes),length(level));
levelGenes(yes,:)=genesrnainfovalue(row(yes),level);
levelGenes(not(yes),:)= mean(levelGenes(yes,:)); % for no data assume the mean of all - this will eliminate the gene from the analysis

% hypnozoites
% [yes,row]=ismember(genes,genesrnainfo);
% sum(yes)
% allval = genesrnainfovalue(row(yes));
% allval(allval>500) = 500;
% hist(allval)
% %
% selectThr = 150;
% levelGenes=nan(length(model.genes),length(level));
% levelGenes(yes,:)=allval;
% indLow = find(levelGenes<selectThr);
% indHigh = find(levelGenes>selectThr);



%check problematic gene
model.genes(not(yes))
% [yes2,row2]=ismember('PBANKA_1118100',genesrnainfo);
% [yes2,row2]=ismember('PBANKA_1118100.1',genesrnainfo);

%%
[yes,row]=ismember(genes,genesrnainfo);
levelGenes_all=zeros(length(model.genes),size(genesrnainfovalue,2));
for i=1:size(genesrnainfovalue,2)
    levelGenes_all(yes,i)=genesrnainfovalue(row(yes),i);
    levelGenes_all(not(yes),i)=mean(levelGenes_all(:,i)); % for no data assume the mean of all - this will eliminate the gene from the analysis
end
