        visualisation=true;
        load 3D_H_store_doppler5_20000;
        
%         if visualisation
%       %visualization of channel estimation from first antenna    
%        figure;
%        subplot(1,2,1);
%        imagesc(abs(H_hard_wnd(:,1:14,1)));
%        xlabel('OFDM Symbol');
%        ylabel('Subcarrier');
%        title('Hard Window Estimate Magnitude');
%        subplot(1,2,2);
%        imagesc(abs(H_idl(:,1:14,1)));
%        xlabel('TTIs');
%        ylabel('Subcarrier');
%        title('Ideal Estimate Magnitude');
%        end

RBs=1;
TX_ant=16;
numUnits = 9; % number of RES block  Conv3D - BN - ReLu
netWidth = 64; % number of convolutional filters
 miniBatchSize=32;
     maxEpochs=50;
%gpuDevice(2);

    
 for i=1:size(H_idl,2)/14
    H_ideal_for_train_cmplx(:,1:14,:,i)=H_idl(:,(i-1)*14+1:i*14,:);
    H_hard_wnd_for_train_cmplx(:,1:14,:,i)=H_hard_wnd(:,(i-1)*14+1:i*14,:);  
 end
%size(H_ideal_for_train_cmplx)
%size(H_hard_wnd_for_train_cmplx)


 for i=1:size(H_idl,2)/14
    H_ideal_for_train_real(:,1:28,1:TX_ant,1,i)=  [real(H_ideal_for_train_cmplx(:,:,1:TX_ant,i)) imag(H_ideal_for_train_cmplx(:,:,1:TX_ant,i))] ;
    H_hard_wnd_train_real(:,1:28,1:TX_ant,1,i)=[real( H_hard_wnd_for_train_cmplx(:,:,1:TX_ant,i)) imag( H_hard_wnd_for_train_cmplx(:,:,1:TX_ant,i))] ; 
 end   
    
 
%  %%%ResNET 3d Architecture without skip connection
%  
%  
% layers = [ ...
%     image3dInputLayer([RBs*12 28 TX_ant 1])
%     
% 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%      batchNormalizationLayer
%     swishLayer 
%     
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     batchNormalizationLayer
%     swishLayer 
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     batchNormalizationLayer
%     reluLayer 
%     convolution3dLayer(3 ,TX_ant, 'Padding', [1,1,1])
%     
%     
%     
%     
%     convolution3dLayer(3 ,1, 'Padding', [1,1,1])
%     regressionLayer];
%     analyzeNetwork(layers);
%  


layers = ResNet3Dgraph(RBs,TX_ant,netWidth,numUnits,"standard");
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
plot(layers)
analyzeNetwork(layers);




% options = trainingOptions('rmsprop', ...
%                          'InitialLearnRate', 3e-4, ...
%                          'SquaredGradientDecayFactor', 0.99,...
%                          'MaxEpochs', 200,...
%                          'MiniBatchSize',128, ...
%                          'Plots', 'training-progress',...
%                          'Verbose',1);


%      maxEpochs=50;
%      l1reg = 0.0001;
%      options = trainingOptions('adam', ...
%         'LearnRateSchedule','piecewise',...
%         'InitialLearnRate', 5e-3, ...
%         'LearnRateDropPeriod',25, ...        
%         'LearnRateDropFactor',0.98, ...
%         'MiniBatchSize',128, ...
%         'MaxEpochs', maxEpochs,...
%         'Plots', 'training-progress',...
%          'Verbose',1, 'VerboseFrequency',1);
     
     
  
     %l1reg = 0.001;
     validationFrequency = floor(2000/miniBatchSize);
     options = trainingOptions('adam', ...        
        'MiniBatchSize',miniBatchSize, ...
        'MaxEpochs', maxEpochs,...
        'Plots', 'training-progress',...
         'Verbose',1, 'VerboseFrequency',1);   
     
     
%      maxEpochs = 100;
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
     
 


        %'L2Regularization',0.0001, ...
        %'L2Regularization',l2reg, ...        
        %'ExecutionEnvironment', 'multi-gpu',...    
        %'Verbose',1,...
        %'VerboseFrequency',10);
 
      %[H_ls_for_train_real_Z,H_ls_for_train_real_mu,H_ls_for_train_real_sigma] = zscore(H_ls_for_train_real);
      %[H_ideal_for_train_real_Z,H_ideal_for_train_real_mu,H_ideal_for_train_real_sigma] = zscore(H_ideal_for_train_real);
      
       Y_residual=H_hard_wnd_train_real-H_ideal_for_train_real;
      
       
       modelDateTime = datestr(now, 'dd-mmm-yyyy-HH-MM-SS');
       %[net, info] = trainNetwork(X_est, Y_residual, layers, options);
       %[net, info] = trainNetwork(H_ls_for_train_real_Z, Y_residual_Z, layers, options);
       [net, info] = trainNetwork(H_hard_wnd_train_real, Y_residual, layers, options);
       save(['trained3dDnCNN-' modelDateTime '-Epoch-' num2str(maxEpochs) '.mat'],'net', 'options');

