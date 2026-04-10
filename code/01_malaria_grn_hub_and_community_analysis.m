%% 01_malaria_grn_hub_and_community_analysis
% This script constructs a correlation-based gene regulatory network (GRN)
% from malaria gene expression data, identifies central hub genes, and
% detects gene communities using hierarchical clustering.
%
% Workflow:
% 1. Load malaria gene expression data
% 2. Filter genes based on nonzero expression
% 3. Normalize expression values using CP10K + log1p
% 4. Filter abnormal cells
% 5. Select highly variable genes
% 6. Build a correlation-based GRN
% 7. Identify top hub genes using degree centrality
% 8. Detect communities using hierarchical clustering
% 9. Save figures and summary tables
%
% Expected input:
%   data/pf_ss2_exp.csv
%
% Outputs:
%   results/Top_Hub_Genes.csv
%   results/Gene_Regulatory_Network.png
%   results/GRN_Dendrogram.png
%   results/Gene_Regulatory_Network_Communities.png

clc;
clear;
close all;

%% ===================== USER SETTINGS =====================
inputFile = fullfile('data', 'pf_ss2_exp.csv');
resultsDir = 'results';

minNonzeroCells = 60;       % Minimum number of nonzero cells required per gene
topExpressedCount = 5;      % Number of top nonzero values used for mean filtering
maxHvGenes = 2000;          % Maximum number of highly variable genes
corrThreshold = 0.7;        % Correlation threshold for GRN construction
topHubCount = 10;           % Number of top hub genes to report
numClusters = 3;            % Number of communities for hierarchical clustering

%% ===================== SETUP =====================
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% ===================== STEP 1: LOAD DATA =====================
dataTable = readtable(inputFile);

% Extract expression matrix and gene names
geneExpression = dataTable{:, 2:end};   % genes x cells
geneNames = string(dataTable{:, 1});

%% ===================== STEP 2: GENE FILTERING =====================
selectedGeneIndices = [];
meanTopExpression = [];

for i = 1:size(geneExpression, 1)
    nonzeroIdx = find(geneExpression(i, :) ~= 0);

    if numel(nonzeroIdx) > minNonzeroCells
        sortedNonzero = sort(geneExpression(i, nonzeroIdx), 'descend');
        selectedGeneIndices(end+1) = i; %#ok<AGROW>
        meanTopExpression(end+1) = mean(sortedNonzero(1:min(topExpressedCount, numel(sortedNonzero)))); %#ok<AGROW>
    end
end

% Sort retained genes by mean of top expressed values
[~, sortIdx] = sort(meanTopExpression, 'descend');
selectedGeneIndices = selectedGeneIndices(sortIdx);

filteredGenes = geneExpression(selectedGeneIndices, :);
filteredGeneNames = geneNames(selectedGeneIndices);

fprintf('Number of genes retained after expression filtering: %d\n', size(filteredGenes, 1));

%% ===================== STEP 3: NORMALIZATION =====================
% Compute total counts per cell
totalCountsPerCell = sum(filteredGenes, 1);
totalCountsPerCell(totalCountsPerCell == 0) = 1;

% Normalize to counts per 10,000
countsPer10k = filteredGenes ./ totalCountsPerCell * 1e4;

% Apply log1p transformation
normalizedData = log1p(countsPer10k);

%% ===================== STEP 4: FILTER ABNORMAL CELLS =====================
cellTotals = sum(filteredGenes, 1);
meanCellTotals = mean(cellTotals);
stdCellTotals = std(cellTotals);

lowerBound = meanCellTotals - stdCellTotals;
upperBound = meanCellTotals + stdCellTotals;

acceptableCells = (cellTotals > lowerBound) & (cellTotals < upperBound);
finalData = normalizedData(:, acceptableCells);

fprintf('Number of cells retained after filtering: %d\n', size(finalData, 2));

%% ===================== STEP 5: SELECT HIGHLY VARIABLE GENES =====================
geneVariance = var(finalData, 0, 2);
numAvailableGenes = numel(geneVariance);

fprintf('Number of genes after cell filtering: %d\n', numAvailableGenes);

numHvGenes = min(maxHvGenes, numAvailableGenes);
[~, varIdx] = sort(geneVariance, 'descend');

hvGenes = finalData(varIdx(1:numHvGenes), :);
hvGeneNames = filteredGeneNames(varIdx(1:numHvGenes));

fprintf('Number of highly variable genes selected: %d\n', size(hvGenes, 1));

%% ===================== STEP 6: PCA =====================
[coeff, score, ~, ~, explained] = pca(hvGenes');

cumExplained = cumsum(explained);
numComponents = find(cumExplained >= 90, 1);

reducedData = score(:, 1:numComponents); %#ok<NASGU>

fprintf('Number of principal components selected: %d\n', numComponents);

%% ===================== STEP 7: BUILD GENE REGULATORY NETWORK =====================
corrMatrix = corr(hvGenes');

adjMatrix = abs(corrMatrix) > corrThreshold;
adjMatrix = adjMatrix - diag(diag(adjMatrix));  % remove self-loops

geneNetwork = graph(adjMatrix, cellstr(hvGeneNames));

% Save GRN figure
fig1 = figure('Visible', 'off');
plot(geneNetwork, 'Layout', 'force');
title('Gene Regulatory Network');
saveas(fig1, fullfile(resultsDir, 'Gene_Regulatory_Network.png'));
close(fig1);

%% ===================== STEP 8: IDENTIFY HUB GENES =====================
degreeCentrality = centrality(geneNetwork, 'degree');

[sortedDegree, idxDegree] = sort(degreeCentrality, 'descend');
topGeneNames = string(geneNetwork.Nodes.Name(idxDegree(1:min(topHubCount, numel(idxDegree)))));
topGeneDegrees = sortedDegree(1:min(topHubCount, numel(sortedDegree)));

hubTable = table(topGeneNames, topGeneDegrees, ...
    'VariableNames', {'Gene', 'DegreeCentrality'});

writetable(hubTable, fullfile(resultsDir, 'Top_Hub_Genes.csv'));

fprintf('Top %d central genes (hubs):\n', min(topHubCount, numel(topGeneNames)));
disp(topGeneNames);

%% ===================== STEP 9: COMMUNITY DETECTION =====================
distMatrix = distances(geneNetwork);

% Replace disconnected distances with a large finite value
finiteDistances = distMatrix(~isinf(distMatrix));
if isempty(finiteDistances)
    error('The network is fully disconnected. Community detection cannot proceed.');
end

maxFiniteDist = max(finiteDistances);
distMatrix(isinf(distMatrix)) = maxFiniteDist + 1;

Y = squareform(distMatrix);
Z = linkage(Y, 'average');

% Save dendrogram
fig2 = figure('Visible', 'off');
dendrogram(Z);
title('Hierarchical Clustering Dendrogram');
saveas(fig2, fullfile(resultsDir, 'GRN_Dendrogram.png'));
close(fig2);

communityIndices = cluster(Z, 'maxclust', numClusters);
geneNetwork.Nodes.Community = communityIndices;

%% ===================== STEP 10: SAVE COMMUNITY NETWORK FIGURE =====================
fig3 = figure('Visible', 'off');
p = plot(geneNetwork, 'Layout', 'force');

p.NodeLabel = geneNetwork.Nodes.Name;
p.MarkerSize = 5;
p.NodeCData = geneNetwork.Nodes.Community;

colormap(jet(numClusters));
colorbar;
title('Gene Regulatory Network with Communities');

saveas(fig3, fullfile(resultsDir, 'Gene_Regulatory_Network_Communities.png'));
close(fig3);

disp('Malaria GRN, hub gene, and community analysis completed successfully.');
