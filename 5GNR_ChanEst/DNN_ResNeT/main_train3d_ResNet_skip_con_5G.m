        visualisation=true;



RBs=20;
TX_ant=4;
numUnits = 15; % number of RES block  Conv3D - BN - ReLu
netWidth = 4*TX_ant; % number of convolutional filters


  

layers = ResNet3Dgraph(RBs,TX_ant,netWidth,numUnits,"standard");
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
plot(layers)
analyzeNetwork(layers);


