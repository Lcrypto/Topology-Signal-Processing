function TDA


rng('default')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%use EPA channel for TDA
%load H_store_EPA.mat;

%use ETU channel for TDA
load H_store_ETU.mat;

  
RBs=1;
TX_ant=8;
 for i=1:size(H_idl,2)/TX_ant
    H_ideal_for_train_cmplx(:,1:TX_ant,1,i)=H_idl(:,(i-1)*TX_ant+1:i*TX_ant);
    H_ls_for_train_cmplx(:,1:TX_ant,1,i)=H_ls(:,(i-1)*TX_ant+1:i*TX_ant);    
 end

 for i=1:size(H_idl,2)/TX_ant
    H_ideal_for_train_real(:,1:TX_ant*2,1,i)=  [real(H_ideal_for_train_cmplx(:,1:TX_ant,1,i)) imag(H_ideal_for_train_cmplx(:,1:TX_ant,1,i))] ;
    H_ls_for_train_real(:,1:TX_ant*2,1,i)=[real( H_ls_for_train_cmplx(:,1:TX_ant,1,i)) imag( H_ls_for_train_cmplx(:,1:TX_ant,1,i))] ; 
 end  
 
%use residual channel 
Y_residual=H_ls_for_train_real-H_ideal_for_train_real;

%use pure channel
%Y_residual=H_ideal_for_train_real;  



B = squeeze(Y_residual(1:RBs*12,1:TX_ant*2,1,1));
data=B;

% B2 = squeeze(Y_residual(1:RBs*12,1:TX_ant*2,1,1:14));
% B3=cat(2,B2(:,:,1),B2(:,:,2),B2(:,:,3),B2(:,:,4),B2(:,:,5),B2(:,:,6),B2(:,:,7),B2(:,:,8),B2(:,:,9),B2(:,:,10),B2(:,:,11),B2(:,:,12),B2(:,:,13),B2(:,:,14));
% data =B3;




% To use instead of Wireless channel, correlation matrix of channel just
% place your correlation matrix in corr_Mat and uncomment lines below

% corr_Mat= 
% [1     0.06  0.23  0.01  0.89
%  0.06  1     0.74  0.01  0.61
%  0.23  0.74  1     0.72  0.03
%  0.01  0.01  0.72  1     0.7 
%  0.89  0.61  0.03  0.7   1   
% ];
% 
% threshold = 0;
%   for i =1: size(corr_Mat,2)
%      for j =1: size(corr_Mat,2)
%       corr_Mat(i,j) = 1 - corr_Mat(i,j);
%         if corr_Mat(i,j) > threshold    
%          threshold = corr_Mat(i,j);
%         end
%      end
%   end
% data =corr_Mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Paths
addpath('SupportFunctions');
addpath('AnalyzeComplex');
addpath('CreateComplex');
addpath('Visualization');
addpath('VRComplex');


%General Parameters

N=size(data,1);
k=1;            %Order of simplices

%Witness Complex Parameters 
%NOTE: Any simplicial complex can be used.  
% The Rips complex is probably the simplest to understand and read about. 
% Zomorodian wrote a great paper about it (see reference above).

R=0.9;          %Connectivity parameters (reduce R to decrease connectivity.)
v=3;            %A witness complex parameter.
n=N;            %nLandmarks (# of points used to aproximate full data set,
                %            this can be either the full set or a subset
                %            chosen using minimax sampling).
                %For large datasets, it makes sense to use a small n for efficiency.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%

data=normalizeData(data);
landmarks=getLandmarksMinMax(data,n,N);

tic
betti=getBettiNumbers(data,landmarks,v,k,R,1);
toc

%Results
for k=1:size(betti,2)
fprintf('Betti%1d: %1d\n',k-1,betti(k));
end

end