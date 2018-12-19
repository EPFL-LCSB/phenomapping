function [expModel, solution, newgrRules, expModelm] = integrateGeneExp(tmodel, ...
    levelGenes, mingrRate, lowPvalue, highPvalue, percLow, percHigh, ...
    numAlt, selectAlt, flagPlot, minmax, geneKO, rxnLevelOUT)
% Integrate gene expression data in the model
%
% USAGE:
%
%       [expModel, solution, newgrRules, expModelm] = integrateGeneExp(tmodel, mingrRate, lowPvalue, highPvalue, percHigh, percLow, numAlt, selectAlt, flagPlot, minmax, geneKO)
%
% INPUTS:
%    tmodel:          model with TFA structure
%    levelGenes:      quantitative gene expression value for all genes in
%                     the model. If lacking data, it is recomended to
%                     assign an average expression value for that genes.
%
% OPTIONAL INPUTS:
%    mingrRate:       min growth to achieve (default = 0.01)
%    lowPvalue:       percentile to assign lowly expressend genes (default
%                     = 20)
%    highPvalue:      percentile to assign highly expressend genes (default
%                     = 80)
%    percLow:         percentage from  (default
%                     = 80)
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    selectAlt:       Indicate what alternative expression profile to
%                     integrate in the model (default = common profile to 
%                     max consistency score (CS))
%    flagPlot:        plot distribution of gene expression (default = 1)
%    minmax:          minmax of net reaction fluxes in two columns 
%                     (default = lb and ub)
%    geneKO:          gene IDs to KO
%
% OUTPUTS:
%    expModel:        model with TFA structure, gene expression 
%                     constraints and information integrated (check field
%                     ExpInfo)
%    solution:        solution structure with information on CS,
%                     alternatives, etc
%    newgrRules:      new GPR rules modified from the original GPR (low 
%                     expressed genes are excluded from complex GPR)
%    expModelm:       model with TFA structure and gene expression 
%                     constraints without a profile integrated
% 
%
% .. Author:
% Daniel F. Hernandez & Vikash Pandey 2015
% Anush Chiappino-Pepe 2017 - integration of expression constraints into 
% the model and organization of the function
%

if (nargin < 3)
    mingrRate = 0.1;
end
if (nargin < 4)
    lowPvalue = 20;
end
if (nargin < 5)
    highPvalue = 80;
end
if (nargin < 6)
    percLow = 2E-5;
end
if (nargin < 7)
    percHigh = 2E-3;
end
if (nargin < 8)
    numAlt = 5;
end
if (nargin < 9)
    selectAlt = 0;
end
if (nargin < 10)
    flagPlot = 1;
end
if (nargin < 11)
    minmax = [];
end
if (nargin < 12)
    geneKO = [];
end
if (nargin < 13)
    rxnLevelOUT = {};
end

if isempty(minmax)
    minmax=[tmodel.lb tmodel.ub];
end

path_save = '/Users/Anush/Documents/gei_altsol.mat';

% Prepare model for gene expression data integration
if flagPlot
    % change gene to rxn expression
%     levels = geneToreaction_levels(tmodel, tmodel.genes, levelGenes, @min, @max );
    lowPer1 = prctile(log(levelGenes),lowPvalue);
    higPer1 = prctile(log(levelGenes),highPvalue);
    % plot distribution
    hist(log(levelGenes),10)
    hold on
    plot([lowPer1; lowPer1], [300*ones(1,length(lowPer1)); ...
        zeros(1,length(lowPer1))])
    plot([higPer1; higPer1], [300*ones(1,length(higPer1)); ...
        zeros(1,length(higPer1))])
    
    ylabel('Number of genes','FontWeight','bold','FontSize',12);
    xlabel('LOG(expression level)','FontWeight','bold','FontSize',12);
    set(gca,'FontSize',10,'FontName','Arial','FontWeight','normal')
    set(gcf, 'Units', 'centimeters','Position', [15 15 10 7.5])
end

% blood
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/14_pbe/data/Pbe_blood_RNAseq.mat')
% load('/Users/Anush/SWITCHdrive/EPFL_Anush/iPfav2/data/2010_Otto_SI4_pfa_data.mat')
