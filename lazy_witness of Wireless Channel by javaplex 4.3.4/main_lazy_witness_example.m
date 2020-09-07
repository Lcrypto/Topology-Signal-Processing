% This script calculates the intervals for a lazy witness complex (Flag
% complex)

clc; clear; close all;
import edu.stanford.math.plex4.*;



%load H_store_EPA.mat;
load H_store_ETU.mat;
 %residual=H_ls-H_idl; % Residual space topology
% data=imag(residual(1:12,1:8));

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
Y_residual=H_ls_for_train_real-H_ideal_for_train_real;
B = squeeze(Y_residual(1:RBs*12,1:TX_ant*2,1,1));
data =B;



max_dimension = 3;
nu = 1;
num_divisions = 10000;


frequency_domain=true;
% each frequency(subcariers domain) is center of ball, use this persistance homology 
% analysis when correlation not exist in space (antennas) like EPA/ETU etc
if frequency_domain
data = B;
num_landmark_points = RBs*12-1; 
end

% each antennas(space domain) is center of ball, use this persistance homology 
% analysis when correlation exist in space (antennas) like EPA/ETU etc
if ~frequency_domain
data = B';
num_landmark_points = size(data,1)-1; 
end

point_cloud=data;
% create a sequential maxmin landmark selector
landmark_selector = api.Plex4.createMaxMinSelector(point_cloud, num_landmark_points);
R = landmark_selector.getMaxDistanceFromPointsToLandmarks()
max_filtration_value = 2 * R;

% create a lazy witness stream
stream = streams.impl.LazyWitnessStream(landmark_selector.getUnderlyingMetricSpace(), landmark_selector, max_dimension, max_filtration_value, nu, num_divisions);
stream.finalizeStream()

% print out the size of the stream
num_simplices = stream.getSize()

% get persistence algorithm over Z/2Z
persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);

% compute the intervals
intervals = persistence.computeIntervals(stream);

% create the barcode plots
options.filename = 'lazySphere';
options.max_filtration_value = max_filtration_value;
options.max_dimension = max_dimension - 1;
plot_barcodes(intervals, options);
