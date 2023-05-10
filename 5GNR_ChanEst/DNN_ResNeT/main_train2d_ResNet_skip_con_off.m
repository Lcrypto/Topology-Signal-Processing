



RBs=1;
TX_ant=8;


layers = [ ...
    imageInputLayer([RBs*12 TX_ant*2 1])
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    reluLayer
    
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
 
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer     
    
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer

    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer    

        convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer  
    convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
    batchNormalizationLayer
    reluLayer 



    
    convolution2dLayer([3 3],1, 'Padding', [1,1])
    regressionLayer];
  

        
      analyzeNetwork(layers);
    
        
