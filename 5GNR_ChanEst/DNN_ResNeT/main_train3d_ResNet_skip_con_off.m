        visualisation=true;



RBs=1;
TX_ant=8;
numUnits = 15; % number of RES block  Conv3D - BN - ReLu
netWidth = 4*TX_ant; % number of convolutional filters


  

layers = [ ...
    image3dInputLayer([RBs*12 28 TX_ant 1])
    
    
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 

    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 

    
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 

    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 

    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 
    convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    batchNormalizationLayer
    reluLayer 


%improve train    
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     
%      batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer
%                                                                                                                                         
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
    
%         convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     batchNormalizationLayer
%     reluLayer 
%     
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     batchNormalizationLayer
%     reluLayer  
%     
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     batchNormalizationLayer
%     reluLayer  
%    
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     batchNormalizationLayer
%     reluLayer  
%     
%     convolution3dLayer(3 ,netWidth, 'Padding', [1,1,1])
%     batchNormalizationLayer
%     reluLayer  
    
    
    convolution3dLayer(3 ,1, 'Padding', [1,1,1])
    regressionLayer];
    analyzeNetwork(layers);




