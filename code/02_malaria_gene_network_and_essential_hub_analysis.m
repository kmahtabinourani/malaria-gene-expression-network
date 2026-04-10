%% 02_malaria_gene_network_and_essential_hub_analysis
% This script analyzes malaria-related gene expression networks and highlights
% essential target genes together with their directly connected neighboring genes.
%
% Workflow:
% 1. Read PDB files from the input folder
% 2. Process protein structures and compute chain-level interaction adjacency
% 3. Read malaria gene expression data
% 4. Extract a target gene name from each PDB filename
% 5. Build a correlation-based gene subnetwork around the target gene
% 6. Highlight the essential gene and its neighboring genes
% 7. Save adjacency matrices, neighboring gene tables, and network figures
%
% Expected inputs:
%   - pdb_inputs/*.pdb
%   - data/pf_ss2_exp.csv
%
% Outputs:
%   - results/<TargetGene>_PPI_AdjacencyMatrix.csv
%   - results/<TargetGene>_NeighboringGenes.csv
%   - results/<TargetGene>_GeneNetwork.png

clc;
clear;
close all;

%% ===================== USER SETTINGS =====================
pdbInputDir = 'pdb_inputs';
expressionFile = fullfile('data', 'pf_ss2_exp.csv');
resultsDir = 'results';

interactionThreshold = 5;   % Angstrom threshold for chain interaction
corrThreshold = 0.8;        % Correlation threshold for gene network edges

%% ===================== SETUP =====================
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

pdbFiles = dir(fullfile(pdbInputDir, '*.pdb'));

if isempty(pdbFiles)
    error('No PDB files were found in the folder: %s', pdbInputDir);
end

%% ===================== LOAD GENE EXPRESSION DATA =====================
dataTable = readtable(expressionFile);
geneExpression = dataTable{:, 2:end};
geneNames = string(dataTable{:, 1});

%% ===================== PROCESS EACH PDB FILE =====================
for kk = 1:length(pdbFiles)

    currentFileName = pdbFiles(kk).name;
    currentFilePath = fullfile(pdbInputDir, currentFileName);

    fprintf('\nProcessing file: %s\n', currentFileName);

    %% -------- Read PDB structure --------
    pdbStruct = pdbread(currentFilePath);

    %% -------- Process protein structure --------
    activeRegions = processProteinStructure(pdbStruct);

    % Extract atomic information
    atoms = pdbStruct.Model.Atom;
    chains = {atoms.chainID};
    uniqueChains = unique(chains);

    positions = [[atoms.X]' [atoms.Y]' [atoms.Z]'];
    atomSeqNumbers = [atoms.resSeq]';

    %#ok<NASGU>
    activeIndices = ismember(atomSeqNumbers, activeRegions);

    %% -------- Build protein-protein interaction adjacency matrix --------
    numChains = numel(uniqueChains);
    ppiAdjMatrix = zeros(numChains);

    for i = 1:numChains
        for j = i+1:numChains

            chainAIndices = strcmp(chains, uniqueChains{i});
            chainBIndices = strcmp(chains, uniqueChains{j});

            chainAPos = positions(chainAIndices, :);
            chainBPos = positions(chainBIndices, :);

            pairwiseDistances = pdist2(chainAPos, chainBPos);

            if any(pairwiseDistances(:) < interactionThreshold)
                ppiAdjMatrix(i, j) = 1;
                ppiAdjMatrix(j, i) = 1;
            end
        end
    end

    %% -------- Save PPI adjacency matrix --------
    [~, baseName, ~] = fileparts(currentFileName);
    writematrix(ppiAdjMatrix, fullfile(resultsDir, [baseName '_PPI_AdjacencyMatrix.csv']));

    %% -------- Extract target gene from file name --------
    targetGene = extractTargetGeneFromFilename(baseName);

    if strlength(targetGene) == 0
        warning('Could not extract target gene from file name: %s. Skipping.', currentFileName);
        continue;
    end

    geneIndex = find(strcmp(geneNames, targetGene), 1);

    if isempty(geneIndex)
        warning('Target gene "%s" was not found in the gene expression table. Skipping.', targetGene);
        continue;
    end

    %% -------- Build correlation-based gene subnetwork --------
    fullCorrMatrix = corr(geneExpression');

    allCorrsToTarget = fullCorrMatrix(geneIndex, :);
    targetEdges = find(abs(allCorrsToTarget) > corrThreshold);

    subNetworkGenes = unique([geneIndex; targetEdges(:)]);
    subAdjMatrix = abs(fullCorrMatrix(subNetworkGenes, subNetworkGenes)) > corrThreshold;

    subGeneNames = geneNames(subNetworkGenes);
    Ggenes = graph(subAdjMatrix, cellstr(subGeneNames));

    %% -------- Identify essential gene and neighboring genes --------
    essentialGeneIndex = find(strcmp(subGeneNames, targetGene), 1);

    nodeColors = repmat([0 0.4470 0.7410], length(subNetworkGenes), 1); % Blue default

    neighboringGenes = strings(0,1);

    if ~isempty(essentialGeneIndex)
        % Highlight essential gene in red
        nodeColors(essentialGeneIndex, :) = [1 0 0];

        % Find directly connected neighboring genes
        connectedEssentialNodes = find(subAdjMatrix(essentialGeneIndex, :) > 0);

        if ~isempty(connectedEssentialNodes)
            % Highlight neighboring genes in yellow
            nodeColors(connectedEssentialNodes, :) = repmat([1 1 0], length(connectedEssentialNodes), 1);

            % Extract neighboring gene names
            neighboringGenes = subGeneNames(connectedEssentialNodes);

            fprintf('Essential Gene: %s\n', targetGene);
            fprintf('Number of Neighboring Genes: %d\n', numel(neighboringGenes));
            disp('Neighboring Genes:');
            disp(neighboringGenes);
        end
    end

    %% -------- Save neighboring genes table --------
    if ~isempty(neighboringGenes)
        neighborTable = table(neighboringGenes, 'VariableNames', {'NeighboringGene'});
        writetable(neighborTable, fullfile(resultsDir, [targetGene '_NeighboringGenes.csv']));
    end

    %% -------- Plot and save gene subnetwork --------
    fig = figure('Visible', 'off');
    h = plot(Ggenes, ...
        'Layout', 'force', ...
        'NodeColor', nodeColors, ...
        'MarkerSize', 7, ...
        'LineWidth', 1.5);

    title(['Gene-Gene Network: ' targetGene ' highlighted']);

    if ~isempty(essentialGeneIndex)
        highlight(h, essentialGeneIndex, 'MarkerSize', 10, 'NodeColor', 'r');
    end

    saveas(fig, fullfile(resultsDir, [targetGene '_GeneNetwork.png']));
    close(fig);

end

disp('Malaria gene network analysis completed successfully.');

%% ===================== LOCAL FUNCTIONS =====================

function activeRegions = processProteinStructure(pdbStruct)
% Identify hydrophobic residues as candidate active regions

    hydrophobicResidues = {'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'MET', 'TRP'};
    atomRecords = pdbStruct.Model.Atom;

    activeRegions = [];
    processedResidues = [];

    for i = 1:length(atomRecords)
        resName = atomRecords(i).resName;
        resSeq = atomRecords(i).resSeq;

        if ismember(resName, hydrophobicResidues) && ~ismember(resSeq, processedResidues)
            activeRegions = [activeRegions; resSeq]; %#ok<AGROW>
            processedResidues = [processedResidues; resSeq]; %#ok<AGROW>
        end
    end

    activeRegions = unique(activeRegions);
end

function targetGene = extractTargetGeneFromFilename(baseName)
% Extract target gene from filename using underscore-separated structure
% Example expected style: something_gene_part_more

    parts = split(string(baseName), '_');

    if numel(parts) >= 3
        % Adjust this rule if your naming pattern differs
        targetGene = strjoin(parts(2:min(3, numel(parts)-1)), '_');

        % If the extracted form still looks unusual, fallback to second token
        if strlength(targetGene) == 0
            targetGene = parts(2);
        end
    elseif numel(parts) >= 2
        targetGene = parts(2);
    else
        targetGene = "";
    end
end
