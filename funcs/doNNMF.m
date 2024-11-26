function [nnmflabels, consistencyScores, H_final] = doNNMF(realdata, clusterRange, numBootstraps, plotTemplates)

% realdata: should be (neuron x time)
% rlabel: neuron x 1 (brain region labels)

% nnmflabels: output labels clustered by NNMF method
% consistencyScores: hight if neurons are assigned to the same cluster across iterations

numNeurons = size(realdata, 1);
numTimePoints = size(realdata, 2);
consistencyScores = zeros(length(clusterRange), 1);

% use parallel computing
if isempty(gcp('nocreate'))
    parpool;
end

% loop cluster num
for idx = 1:length(clusterRange)
    k = clusterRange(idx);
    disp(['Processing cluster num = ' num2str(k)]);

    clusterAssignments = zeros(numNeurons, numBootstraps);

    % loop bootstrap
    parfor b = 1:numBootstraps
        % train-test split 5:5
        randIndices = randperm(numNeurons);
        trainIndices = randIndices(1:round(numNeurons/2));
        testIndices = randIndices(round(numNeurons/2)+1:end);
        trainingData = realdata(trainIndices, :);
        testData = realdata(testIndices, :);

        % run NNMF on the training data
        % W: neuron x k, H: k x time
        [W_train, H_train] = nnmf(trainingData, k);

        % apply to the test data
        % lsqnonneg: linear least squares with nonnegativity constraints
        % X = LSQNONNEG(C,d) returns the vector X that minimizes NORM(d-C*X)
        W_test = zeros(length(testIndices), k);        
        for i = 1:length(testIndices)
            W_test(i, :) = lsqnonneg(H_train', testData(i, :)')';
        end
        % get cluster label based on the maximum weights
        [~, clusterLabels] = max(W_test, [], 2);

        tempClusterAssignments = zeros(numNeurons, 1);
        tempClusterAssignments(testIndices) = clusterLabels;
        clusterAssignments(:, b) = tempClusterAssignments;
    end
    allNeuronClusters = clusterAssignments;

    neuronConsistency = zeros(numNeurons, 1);
    % calculate consistency scores
    parfor n = 1:numNeurons
        % >0 to get iteration where the neuron was tested
        clusters = allNeuronClusters(n, allNeuronClusters(n, :) > 0);
        if ~isempty(clusters)
            % get mode and its proportion
            modeCluster = mode(clusters);
            consistency = sum(clusters == modeCluster) / length(clusters);
            neuronConsistency(n) = consistency;
        end
    end
    consistencyScores(idx) = mean(neuronConsistency);
    disp(['Mean consistency for cluster num = ' num2str(k) ': ' num2str(consistencyScores(idx))]);
end

% decide the best num of clusters
[~, bestKIdx] = max(consistencyScores);
bestK = clusterRange(bestKIdx);
disp(['Optimal number of clusters: ' num2str(bestK)]);

% do NNMF based on the bestK
[W_final, H_final] = nnmf(realdata, bestK);

% get the final results
[~, nnmflabels] = max(W_final, [], 2);

% show templates
if plotTemplates
    plot(H_final');
    title('Template Activity Patterns');
end

end