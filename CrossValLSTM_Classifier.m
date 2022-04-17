% This script contains the code to load data, create a train/test split and
% train the LSTM model
% Author Shane Malone, Student number 17308336

%load('FeatureSets\FusedFeatureSetMIMIC20Beats.mat')

%% 

numFeatures = 1;
numHiddenUnits = 40;
numClasses = 2;
%miniBatchSize = 1024;
miniBatchSize = 4096;

Features = feature';


c = cvpartition(length(Features), 'HoldOut', 0.3);
trainIndices = training(c);
testIndices = test(c);

XTrain = Features(trainIndices);
trainLabels = AF(trainIndices);
XValidation = Features(testIndices);
valLabels = AF(testIndices);

YTrain = cell2mat(trainLabels);
YTrain = categorical(YTrain)';

YValidation = cell2mat(valLabels);
YValidation = categorical(YValidation)';


count = 0;
for i=1:length(YTrain)
    if trainLabels{i} == 1
        count = count +1;
    end
end
AFibRatio = count/length(YTrain);
nonAfibRatio = 1 - AFibRatio;

classWeights = [1/nonAfibRatio, 1/AFibRatio]; % Inverse of count per classes
classes = ["0", "1"];

count = 0;
for i=1:length(YValidation)
    if valLabels{i} == 1
        count = count +1;
    end
end
ValAFibRatio = count/length(YValidation);
ValnonAfibRatio = 1 - ValAFibRatio;
%% 

layers = [ ...
    sequenceInputLayer(numFeatures)
    bilstmLayer(numHiddenUnits,'OutputMode','sequence')
    globalMaxPooling1dLayer
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer('Classes', classes, 'ClassWeights',classWeights)
    %classificationLayer
    ];

options = trainingOptions('adam', ...
    'MaxEpochs',250, ...
    'GradientThreshold',0.5, ...
    'MiniBatchSize',miniBatchSize, ...
    'Verbose',0, ...
    'ValidationData', {XValidation,YValidation}, ...
    'OutputNetwork','best-validation-loss', ...
    'Plots','training-progress', ...
    'ExecutionEnvironment','cpu');

net = trainNetwork(XTrain, YTrain, layers, options);
%% 

yPred = classify(net, XValidation);
%figure
%cm = confusionchart(YValidation,yPred);

tp = sum((yPred == categorical(1)) & (YValidation == categorical(1)));
fp = sum((yPred == categorical(1)) & (YValidation == categorical(0)));
tn = sum((yPred == categorical(0)) & (YValidation == categorical(0)));
fn = sum((yPred == categorical(0)) & (YValidation == categorical(1)));

sensitivity = tp/(tp + fn);
specificity = tn/(tn + fp);
precision = tp / (tp + fp);
FPR = fp/(tn+fp);
Accuracy = (tp+tn)./(tp+fp+tn+fn);
recall = tp / (tp + fn);
F1 = (2 * precision * recall) / (precision + recall);