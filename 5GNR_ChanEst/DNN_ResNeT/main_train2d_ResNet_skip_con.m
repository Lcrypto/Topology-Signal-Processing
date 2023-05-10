% load dumped channel data 
% load 2D_H_store_doppler5_20000.mat;
% RBs=size(H_idl2d,1)/12;
% TX_ant=16;


TX_ant=8;
numUnits = 21; % number of RES block  Conv2D - BN - ReLu
netWidth = 2*TX_ant; % number of convolutional filters

miniBatchSize=32;
     maxEpochs=100;

RBs=1;




layers = ResNetgraph(RBs,TX_ant,netWidth,numUnits,"standard");
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
plot(layers)
        
      analyzeNetwork(layers);
    


     
     
     
     