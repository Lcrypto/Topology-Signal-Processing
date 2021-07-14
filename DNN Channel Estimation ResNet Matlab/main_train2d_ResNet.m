% load dumped channel data 
% load 2D_H_store_doppler5_20000.mat;
% RBs=size(H_idl2d,1)/12;
% TX_ant=16;


load 2D_H_store_doppler5_2000.mat; 
RBs=size(H_idl2d,1)/12;
TX_ant=64;

numUnits = 9; % number of RES block  Conv2D - BN - ReLu
netWidth = 64; % number of convolutional filters

miniBatchSize=32;
     maxEpochs=100;

%gpuDevice(0);
sft_window= false ; %true for train use sft_wnd, false ls or hw avg window depend from platform
H_idl=H_idl2d;
H_ls=H_hard_wnd2d;
if sft_window
 for i=1:size(H_idl_2,2)/TX_ant
    H_ideal_for_train_cmplx(:,1:TX_ant,1,i)=H_idl_2(:,(i-1)*TX_ant+1:i*TX_ant);
    H_ls_for_train_cmplx(:,1:TX_ant,1,i)=H_sft_wnd_l(:,(i-1)*TX_ant+1:i*TX_ant);    
 end

 for i=1:size(H_idl_2,2)/TX_ant
    H_ideal_for_train_real(:,1:TX_ant*2,1,i)=  [real(H_ideal_for_train_cmplx(:,1:TX_ant,1,i)) imag(H_ideal_for_train_cmplx(:,1:TX_ant,1,i))] ;
    H_ls_for_train_real(:,1:TX_ant*2,1,i)=[real( H_ls_for_train_cmplx(:,1:TX_ant,1,i)) imag( H_ls_for_train_cmplx(:,1:TX_ant,1,i))] ; 
 end  
else
    
%  for i=1:size(H_idl,2)/TX_ant
%     H_ideal_for_train_cmplx(:,1:TX_ant,1,i)=H_idl(:,(i-1)*TX_ant+1:i*TX_ant);
%     H_ls_for_train_cmplx(:,1:TX_ant,1,i)=H_ls(:,(i-1)*TX_ant+1:i*TX_ant);    
%  end
 
 H_ideal_for_train_cmplx=zeros(size(H_idl,1),TX_ant,size(H_idl,2)/TX_ant);
 H_ls_for_train_cmplx=zeros(size(H_idl,1),TX_ant,size(H_idl,2)/TX_ant);
  for i=1:size(H_idl,2)/TX_ant
    H_ideal_for_train_cmplx(:,1:TX_ant,i)=H_idl(:,(i-1)*TX_ant+1:i*TX_ant);
    H_ls_for_train_cmplx(:,1:TX_ant,i)=H_ls(:,(i-1)*TX_ant+1:i*TX_ant);    
  end
 
H_ideal_for_train_real=zeros(size(H_ideal_for_train_cmplx,1),TX_ant*2,1,size(H_idl,2)/TX_ant);
H_ls_for_train_real=zeros(size(H_ideal_for_train_cmplx,1),TX_ant*2,1,size(H_idl,2)/TX_ant);
 for i=1:size(H_idl,2)/TX_ant
    H_ideal_for_train_real(:,1:TX_ant*2,1,i)=  [real(H_ideal_for_train_cmplx(:,1:TX_ant,i)) imag(H_ideal_for_train_cmplx(:,1:TX_ant,i))] ;
    H_ls_for_train_real(:,1:TX_ant*2,1,i)=[real( H_ls_for_train_cmplx(:,1:TX_ant,i)) imag( H_ls_for_train_cmplx(:,1:TX_ant,i))] ; 
 end   
    
end

%%%ResNET Architecture without skip connection


%architecture from
%Zhang, K., W. Zuo, Y. Chen, D. Meng, and L. Zhang, "Beyond a Gaussian Denoiser: Residual Learning of Deep CNN for Image Denoising." IEEE® Transactions on Image Processing. Feb 2017.

% layers = [ ...
%     imageInputLayer([RBs*12 TX_ant*2 1])
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     reluLayer('Name','relu_1')  
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_2')  
%     additionLayer(2,'Name','add_1')
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_3') 
%     additionLayer(2,'Name','add_2')
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_3')  
%     additionLayer(2,'Name','add_3')
%     
%      convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_4') 
%     additionLayer(2,'Name','add_4')
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_5') 
%     additionLayer(2,'Name','add_5')    
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_6') 
%     additionLayer(2,'Name','add_6')   
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_7') 
%     additionLayer(2,'Name','add_7')   
% 
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_8') 
%     additionLayer(2,'Name','add_8')      
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_9') 
%     additionLayer(2,'Name','add_9')     
%     
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_9') 
%     additionLayer(2,'Name','add_9')    
%     
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_10') 
%     additionLayer(2,'Name','add_10')      
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_11') 
%     additionLayer(2,'Name','add_11')    
% 
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_12') 
%     additionLayer(2,'Name','add_12')  
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_13') 
%     additionLayer(2,'Name','add_13')      
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_14') 
%     additionLayer(2,'Name','add_14')     
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_15') 
%     additionLayer(2,'Name','add_15')  
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_16') 
%     additionLayer(2,'Name','add_16')  
%     
%     convolution2dLayer([3 3],TX_ant*2, 'Padding', [1,1])
%     batchNormalizationLayer
%     reluLayer('Name','relu_17') 
%     additionLayer(2,'Name','add_17')   
%     
%     
%     convolution2dLayer([3 3],1, 'Padding', [1,1])
%     regressionLayer];
  





layers = ResNetgraph(RBs,TX_ant,netWidth,numUnits,"standard");
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
plot(layers)
        
      analyzeNetwork(layers);
    
        



%          maxEpochs=5;
%          initLearningRate=0.05;
%          l2reg = 0.0001;
%          batchSize=64;
%          options = trainingOptions('sgdm', ...
%         'Momentum',0.9,...
%         'InitialLearnRate', initLearningRate, ...
%         'GradientThresholdMethod','absolute-value', ...        
%         'GradientThreshold',0.005, ...
%         'L2Regularization',l2reg, ...
%         'MiniBatchSize',batchSize, ...
%         'MaxEpochs', maxEpochs,...
%         'Plots', 'training-progress',...
%          'Verbose',1,'VerboseFrequency',1);
     
     
     
  %          'L2Regularization',l2reg, ...        
%        'ExecutionEnvironment', 'multi-gpu',...    
        %'Verbose',1,...
        %'VerboseFrequency',10);
        
%         maxEpochs = 100;
% epochIntervals = 1;
% initLearningRate = 0.1;
% learningRateFactor = 0.1;
% l2reg = 0.0001;
% miniBatchSize = 64;
% options = trainingOptions('sgdm', ...
%     'Momentum',0.9, ...
%     'InitialLearnRate',initLearningRate, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropPeriod',10, ...
%     'LearnRateDropFactor',learningRateFactor, ...
%     'L2Regularization',l2reg, ...
%     'MaxEpochs',maxEpochs, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'GradientThresholdMethod','l2norm', ...
%     'GradientThreshold',0.01, ...
%     'Plots','training-progress', ...
%     'Verbose',false);
        
     
     
     %l1reg = 0.001;

     options = trainingOptions('adam', ...        
        'MiniBatchSize',miniBatchSize, ...
        'MaxEpochs', maxEpochs,...
        'Plots', 'training-progress',...
         'Verbose',1, 'VerboseFrequency',10);
     
     
     
     
     
        
      [H_ls_for_train_real_Z,H_ls_for_train_real_mu,H_ls_for_train_real_sigma] = zscore(H_ls_for_train_real);
      [H_ideal_for_train_real_Z,H_ideal_for_train_real_mu,H_ideal_for_train_real_sigma] = zscore(H_ideal_for_train_real);
      
       Y_residual=H_ls_for_train_real-H_ideal_for_train_real;
       Y_residual_Z=H_ls_for_train_real_Z-H_ideal_for_train_real_Z;
       
       modelDateTime = datestr(now, 'dd-mmm-yyyy-HH-MM-SS');
       %[net, info] = trainNetwork(X_est, Y_residual, layers, options);
       %[net, info] = trainNetwork(H_ls_for_train_real_Z, Y_residual_Z, layers, options);
       [net, info] = trainNetwork(H_ls_for_train_real, Y_residual, layers, options);
       save(['trainedDnCNN-' modelDateTime '-Epoch-' num2str(maxEpochs) '.mat'],'net', 'options');

