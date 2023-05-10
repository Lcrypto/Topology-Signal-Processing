

function lgraph = ResNet3Dgraph(RBs, TX_ant, netWidth,numUnits,unitType)


% Check inputs
assert(numUnits > 0 && mod(numUnits,3) == 0 ,...
    "Number of convolutional units must be an integer multiple of 3.");
unitsPerStage = numUnits/3;

if unitType == "standard"
    convolutionalUnit = @standardConvolutionalUnit;
elseif unitType == "bottleneck"
    convolutionalUnit = @bottleneckConvolutionalUnit;
else
    error("Residual block type must be either ""standard"" or ""bottleneck"".")
end


%% Create Main Network Branch

%
% Input section. Add the input layer and the first convolutional layer.
layers = [
    image3dInputLayer([RBs*12 28 TX_ant 1],'Name','input')
    %image3dInputLayer([RBs*12 TX_ant 28 1],'Name','input')
    convolution3dLayer(3,netWidth,'Padding','same','Name','convInp')
    batchNormalizationLayer('Name','BNInp')
    reluLayer('Name','reluInp')];

% Stage one. Activation size is 32-by-32.
for i = 1:unitsPerStage
    layers = [layers
        convolutionalUnit(netWidth,['S1U' num2str(i) '_'])
        additionLayer(2,'Name',['add1' num2str(i)])
        reluLayer('Name',['relu1' num2str(i)])];
end


for i = 1:unitsPerStage
   
    layers = [layers
        convolutionalUnit(netWidth,['S2U' num2str(i) '_'])
        additionLayer(2,'Name',['add2' num2str(i)])
        reluLayer('Name',['relu2' num2str(i)])];
end


for i = 1:unitsPerStage
  
    layers = [layers
        convolutionalUnit(netWidth,['S3U' num2str(i) '_'])
        additionLayer(2,'Name',['add3' num2str(i)])
        reluLayer('Name',['relu3' num2str(i)])];
end

% Output section.
layers = [layers
         convolution3dLayer(3 ,1,'Padding', [1,1,1],'Name','conv_final')
         regressionLayer('Name','regressOutput')];
    

lgraph = layerGraph(layers);


%% Add shortcut connections
% Add shortcut connection around the convolutional units. Most shortcuts
% are identity connections.
for i = 1:unitsPerStage-1
    lgraph = connectLayers(lgraph,['relu1' num2str(i)],['add1' num2str(i+1) '/in2']);
    lgraph = connectLayers(lgraph,['relu2' num2str(i)],['add2' num2str(i+1) '/in2']);
    lgraph = connectLayers(lgraph,['relu3' num2str(i)],['add3' num2str(i+1) '/in2']);
end

% Shortcut connection from input section to first stage. If unitType equals
% "bottleneck", then the shortcut connection must upsample the channel
% dimension from netWidth to netWidth*4.

    lgraph = connectLayers(lgraph,'reluInp','add11/in2');


    numF =  netWidth;

skip1 = [convolution3dLayer(1,numF,'Name','skipConv1')
    batchNormalizationLayer('Name','skipBN1')];
lgraph = addLayers(lgraph,skip1);
lgraph = connectLayers(lgraph,['relu1' num2str(unitsPerStage)],'skipConv1');
lgraph = connectLayers(lgraph,'skipBN1','add21/in2');


  numF =  netWidth;

skip2 = [convolution3dLayer(1,numF,'Name','skipConv2')
    batchNormalizationLayer('Name','skipBN2')];
lgraph = addLayers(lgraph,skip2);
lgraph = connectLayers(lgraph,['relu2' num2str(unitsPerStage)],'skipConv2');
lgraph = connectLayers(lgraph,'skipBN2','add31/in2');

return


end

%%
% layers = standardConvolutionalUnit(numF,stride,tag) creates a standard
% convolutional unit, containing two 3-by-3 convolutional layers with numF
% filters and a tag for layer name assignment.
function layers = standardConvolutionalUnit(numF,tag)
layers = [
    convolution3dLayer(3,numF,'Padding','same','Name',[tag,'conv1'])
    batchNormalizationLayer('Name',[tag,'BN1'])
    reluLayer('Name',[tag,'relu1'])
    convolution3dLayer(3,numF,'Padding','same','Name',[tag,'conv2'])
    batchNormalizationLayer('Name',[tag,'BN2'])];
end


