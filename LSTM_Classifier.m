load('FeatureSet.mat')

numFeatures = 1;
numHiddenUnits = 100;
numClasses = 2;
miniBatchSize = 10;

XTrain = feature';

%YTrain = cellfun(@categorical, AFperSample, 'uniform', 0)';
YTrain = cell2mat(AF);
YTrain = categorical(YTrain)';


layers = [ ...
    sequenceInputLayer(numFeatures)
    bilstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',80, ...
    'GradientThreshold',7, ...
    'MiniBatchSize',miniBatchSize, ...
    'Verbose',0, ...
    'Plots','training-progress');

net = trainNetwork(XTrain, YTrain, layers, options);