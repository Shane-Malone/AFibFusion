%load('FeatureSetMIMIC20BeatsFixed.mat')
%load('FusedFeatureSetMIMIC20Beats.mat')
load('FeatureSetMIMIC20Beats_Noisy_0dB.mat')
%load('FusedFeatureSetMIMIC20Beats_Noisy_10dB.mat')
%load('FeatureSetMIT20Beats.mat')
%load('PPGFeatureSetMIMIC20Beats.mat')
%% 

numFeatures = 1;
numHiddenUnits = 40;
numClasses = 2;
miniBatchSize = 1024;

Features = feature';
Labels = cell2mat(AF);
%YTrain = cell2mat(AFperSample);
Labels = categorical(Labels)';

%YTrain = AFperSample';


valStart = 1; % Extract valStart->end for validation
valEnd = 94157;
trainStart = 94158;
trainEnd = 285325;

% Split out val data
XValidation = Features(valStart:valEnd);
YValidation = Labels(valStart:valEnd);
XTrain = Features(trainStart:trainEnd);
YTrain = Labels(trainStart:trainEnd);

count = 0;
for i=1:length(YTrain)
    if AF{i+trainStart-1} == 1
        count = count +1;
    end
end
AFibRatio = count/length(YTrain);
nonAfibRatio = 1 - AFibRatio;

classWeights = [1/nonAfibRatio, 1/AFibRatio]; % Inverse of count per classes
classes = ["0", "1"];

count = 0;
for i=1:length(YValidation)
    if AF{i+valStart-1} == 1
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
    ]

options = trainingOptions('adam', ...
    'MaxEpochs',150, ...
    'GradientThreshold',0.5, ...
    'MiniBatchSize',miniBatchSize, ...
    'Verbose',0, ...
    'ValidationData', {XValidation,YValidation}, ...
    'Plots','training-progress', ...
    'ExecutionEnvironment','cpu');

net = trainNetwork(XTrain, YTrain, layers, options);
%% 

yPred = classify(net, XValidation);
figure
cm = confusionchart(YValidation,yPred);