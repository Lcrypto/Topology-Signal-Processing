clear;
%% NR PDSCH Throughput 
% This example shows how to measure the physical downlink shared channel
% (PDSCH) throughput of a 5G New Radio (NR) link, as defined by the 3GPP NR
% standard. The example implements PDSCH and downlink transport channel
% DL-SCH. The transmitter model includes PDSCH demodulation reference
% symbols (DM-RS) and synchronization signal (SS) bursts. The example
% supports both clustered delay line (CDL) and tapped delay line (TDL)
% propagation channels. You can perform perfect or practical
% synchronization and channel estimation. To reduce the total simulation
% time, you can execute the SNR points in the SNR loop in parallel by using
% the Parallel Computing Toolbox(TM).

%% Copyright 2017-2019 The MathWorks, Inc.



%% Introduction
% This example measures the PDSCH throughput of a 5G link, as defined by
% the 3GPP NR standard [ <#14 1> ], [ <#14 2> ], [ <#14 3> ], [ <#14 4> ].
% 
% The following 5G NR features are modeled:
% 
% * DL-SCH transport channel coding 
% * PDSCH and PDSCH DM-RS generation
% * SS burst generation (PSS/SSS/PBCH/PBCH DM-RS)
% * Variable subcarrier spacing and frame numerologies (2^n * 15kHz) for
% normal and extended cyclic prefix
% * TDL and CDL propagation channel models
% 
% Other features of the simulation are:
% 
% * PDSCH precoding using SVD
% * CP-OFDM modulation
% * Slot wise and non slot wise PDSCH and DM-RS mapping
% * SS burst generation (cases A-E, SS/PBCH block bitmap control)
% * Perfect or practical synchronization and channel estimation
% * HARQ operation with 16 processes
% * The example uses a single bandwidth part across the whole carrier
% 
% The figure below shows the processing chain implemented. For clarity, the 
% DM-RS and SS burst generation have been omitted. 
% 
% <<../PDSCHLinkExampleProcessingChain.png>>
% 
% This example supports operation with one or two codewords, dependent on
% the number of layers selected for the transmission. A single precoding
% matrix for the whole PDSCH allocation is determined using SVD by
% averaging the channel estimate across all allocated PDSCH PRBs.
% Therefore, for large PDSCH allocations, i.e. occupying a wide bandwidth,
% the single precoding matrix may not be well matched to the channel across
% all frequencies, resulting in performance degradation. There is no
% beamforming on the SS/PBCH blocks in the SS burst.
%
% To reduce the total simulation time, you can use the Parallel Computing
% Toolbox to execute the SNR points of the SNR loop in parallel.

%% Simulation Length and SNR Points
% Set the length of the simulation in terms of the number of 10ms frames.
% A large number of NFrames should be used to produce meaningful throughput
% results. Set the SNR points to simulate. The SNR for each layer is
% defined per RE, and it includes the effect of signal and noise across all
% antennas.
stream = RandStream('mt19937ar','seed',2);
RandStream.setGlobalStream(stream);
randn

%% Channel Estimator Configuration
% The logical variable |perfectChannelEstimator| controls channel
% estimation and synchronization behavior. When set to |true|, perfect
% channel estimation and synchronization is used. Otherwise, practical
% channel estimation and synchronization is used, based on the values of
% the received PDSCH DM-RS.



%For perfect channel estimation use settings

%perfectChannelEstimator = true;
%inference=0; 
%inference2d=0; 
%collect_data=0; 
%collect_data2d=0;

%For 3D DNN

%perfectChannelEstimator = false;
%inference=1; 
%inference2d=0; 
%collect_data=0; 
%collect_data2d=0;

%For 2D DNN

%perfectChannelEstimator = false;
%inference=0; 
%inference2d=1; 
%collect_data=0; 
%collect_data2d=0;

%For collecting data for Train DNNs 

%perfectChannelEstimator = false;
%inference=0; 
%inference2d=0; 
%collect_data=1; 
%collect_data2d=1;



perfectChannelEstimator = false;
inference=0; % use  3D convolutional Neural Network channel estimation from Hard Window
inference2d=0; % use  2D convolutional Neural Network channel estimation from Hard Window
%use inference=1 only when perfectChannelEstimator = false
perfect_noise_estimation=0; %use perfect noise power for scale, false hard-window estimation
simParameters = [];             % Clear simParameters variable
simParameters.NFrames = 20000;      % Number of frames
%simParameters.SNRIn = [-8 -7 -6 -5 -4 -3 -2 -1 0 ]; % SNR range (dB)
%simParameters.SNRIn = [-2.4]; % SNR range (dB)
simParameters.SNRIn = [-6]; % SNR range (dB)
collect_data=0; %Collect Ideal Channel estimation and Hard_WNd channel estimation for 
%3D Convolutional Network Training
collect_data2d=1; % collect data for 2D Channel Estimation for if  collect_data=1
save_period=10000; %store in file each save_period OFDMs 
global H_idl;
global H_hard_wnd;

H_idl=[];
H_hard_wnd=[];


%% gNodeB and PDSCH Configuration
% Set the key parameters of the simulation. These include:
% 
% * The bandwidth in resource blocks (12 subcarriers per resource block).
% * Subcarrier spacing: 15, 30, 60, 120, 240 (kHz)
% * Cyclic prefix length: normal or extended
% * Cell ID
% * Number of transmit and receive antennas
% 
% A substructure containing the DL-SCH and PDSCH parameters is also
% specified. This includes:
% 
% * Target code rate
% * Allocated resource blocks (PRBSet)
% * Modulation scheme: 'QPSK', '16QAM', '64QAM', '256QAM'
% * Number of layers
% * PDSCH mapping type
% * DM-RS configuration parameters
% 
% Other simulation wide parameters are:
% 
% * Propagation channel model: 'TDL' or 'CDL'
% * SS burst configuration parameters. Note that the SS burst generation
% can be disabled by setting the |SSBTransmitted| field to [0 0 0 0].

% Set waveform type and PDSCH numerology (SCS and CP type)
%simParameters.NRB = 20;                  % Bandwidth in number of resource blocks (51RBs at 30kHz SCS for 20MHz BW)
simParameters.NRB = 20;  
simParameters.SubcarrierSpacing = 30;    % 15, 30, 60, 120, 240 (kHz)
simParameters.CyclicPrefix = 'Normal';   % 'Normal' or 'Extended'
simParameters.NCellID = 1;               % Cell identity

% DL-SCH/PDSCH parameters
simParameters.PDSCH.PRBSet = 0:simParameters.NRB-1; % PDSCH PRB allocation
%simParameters.PDSCH.PRBSet = [0,1]; % PDSCH PRB allocation
%simParameters.PDSCH.PRBSet = [0,1,2,3]; % PDSCH PRB allocation


simParameters.PDSCH.SymbolSet = 0:13;           % PDSCH symbol allocation in each slot
simParameters.PDSCH.EnableHARQ = false;          % Enable/disable HARQ, if disabled, single transmission with RV=0, i.e. no retransmissions

simParameters.PDSCH.NLayers = 1;                % Number of PDSCH layers
simParameters.NTxAnts = 4;                      % Number of PDSCH transmission antennas
if simParameters.PDSCH.NLayers > 4              % Multicodeword transmission
    simParameters.NumCW = 2;                        % Number of codewords
    simParameters.PDSCH.TargetCodeRate = [490 490]./1024; % Code rate used to calculate transport block sizes
    simParameters.PDSCH.Modulation = {'16QAM','16QAM'};   % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.NRxAnts = 4;                            % Number of UE receive antennas
else
    simParameters.NumCW = 1;                        % Number of codewords
    simParameters.PDSCH.TargetCodeRate = 490/1024;  % Code rate used to calculate transport block sizes
    simParameters.PDSCH.Modulation = '16QAM';       % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.NRxAnts = 4;                      % Number of UE receive antennas
end

% DM-RS and antenna port configuration (TS 38.211 Section 7.4.1.1)
simParameters.PDSCH.PortSet = 0:simParameters.PDSCH.NLayers-1; % DM-RS ports to use for the layers
simParameters.PDSCH.PDSCHMappingType = 'A';     % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
simParameters.PDSCH.DMRSTypeAPosition = 2;      % Mapping type A only. First DM-RS symbol position (2,3)
simParameters.PDSCH.DMRSLength = 1;             % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
simParameters.PDSCH.DMRSAdditionalPosition = 0; % Additional DM-RS symbol positions (max range 0...3)
simParameters.PDSCH.DMRSConfigurationType = 2;  % DM-RS configuration type (1,2)
simParameters.PDSCH.NumCDMGroupsWithoutData = 0;% CDM groups without data
simParameters.PDSCH.NIDNSCID = 1;               % Scrambling identity (0...65535)
simParameters.PDSCH.NSCID = 0;                  % Scrambling initialization (0,1)
% Reserved PRB patterns (for CORESETs, forward compatibility etc)
simParameters.PDSCH.Reserved.Symbols = [];      % Reserved PDSCH symbols
simParameters.PDSCH.Reserved.PRB = [];          % Reserved PDSCH PRBs
simParameters.PDSCH.Reserved.Period = [];       % Periodicity of reserved resources
% PDSCH resource block mapping (TS 38.211 Section 7.3.1.6)
simParameters.PDSCH.VRBToPRBInterleaving = 0;   % Disable interleaved resource mapping

% Define the propagation channel type
simParameters.ChannelType = 'CDL'; % 'CDL' or 'TDL'

% SS burst configuration
simParameters.SSBurst.BlockPattern = 'Case B';    % 30kHz subcarrier spacing
%simParameters.SSBurst.SSBTransmitted = [0 1 0 1]; % Bitmap indicating blocks transmitted in the burst
simParameters.SSBurst.SSBTransmitted = [0 0 0 0];
%simParameters.SSBurst.SSBPeriodicity = 20;        % SS burst set periodicity in ms (5, 10, 20, 40, 80, 160)
simParameters.SSBurst.SSBPeriodicity = 160;
validateNLayers(simParameters);

%% 
% Create gNodeB configuration structure 'gnb' and PDSCH configuration
% structure 'pdsch'
gnb = simParameters;
pdsch = simParameters.PDSCH;

% Specify additional required fields for PDSCH
pdsch.RNTI = 1;

Xoh_PDSCH = 0;     % The Xoh-PDSCH overhead value is taken to be 0 here

%% Propagation Channel Model Configuration
% Create the channel model object. Both CDL and TDL channel models are
% supported [ <#14 6> ].

snrIn = simParameters.SNRIn;
nTxAnts = simParameters.NTxAnts;
nRxAnts = simParameters.NRxAnts;
channelType = simParameters.ChannelType;

if strcmpi(channelType,'CDL')
    channel = nrCDLChannel; % CDL channel object
    
    % Use CDL-C model (Urban macrocell model)
    channel.DelayProfile = 'CDL-C';
    channel.DelaySpread = 300e-9;
    channel.MaximumDopplerShift=5;
    channel.CarrierFrequency=3.5e9; %3.5 GHz 
    % Turn the overall number o;;f antennas into a specific antenna panel
    % array geometry
    [channel.TransmitAntennaArray.Size, channel.ReceiveAntennaArray.Size] = ...
        hArrayGeometry(nTxAnts,nRxAnts);
        
elseif strcmpi(channelType,'TDL')
    channel = nrTDLChannel; % TDL channel object
    % Set the channel geometry
    channel.DelayProfile = 'TDL-C';
    channel.DelaySpread = 300e-9;
    channel.CarrierFrequency=3.5e9; % 
    channel.NumTransmitAntennas = nTxAnts;
    channel.NumReceiveAntennas = nRxAnts;
else
    error('ChannelType parameter field must be either CDL or TDL.');
end

%% 
% The sampling rate for the channel model is set using the value returned 
% from <matlab:help('hOFDMInfo') hOFDMInfo>.

waveformInfo = hOFDMInfo(gnb);
channel.SampleRate = waveformInfo.SamplingRate;

%% 
% Get the maximum number of delayed samples by a channel multipath
% component. This is calculated from the channel path with the largest
% delay and the implementation delay of the channel filter. This is
% required later to flush the channel filter to obtain the received signal.

chInfo = info(channel);
maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;

%% Reserve PDSCH Resources Corresponding to SS burst
% This section shows how to reserve resources for the transmission of the 
% SS burst.

% Create SS burst information structure
ssburst = simParameters.SSBurst;
ssburst.NCellID = gnb.NCellID;
ssburst.SampleRate = waveformInfo.SamplingRate;

ssbInfo = hSSBurstInfo(ssburst);

% Map the occupied subcarriers and transmitted symbols of the SS burst
% (defined in the SS burst numerology) to PDSCH PRBs and symbols in the
% data numerology
[mappedPRB,mappedSymbols] = mapNumerology(ssbInfo.OccupiedSubcarriers,ssbInfo.OccupiedSymbols,ssbInfo.NRB,gnb.NRB,ssbInfo.SubcarrierSpacing,gnb.SubcarrierSpacing);

% Configure the PDSCH to reserve these resources so that the PDSCH
% transmission does not overlap the SS burst
reservation.Symbols = mappedSymbols;
reservation.PRB = mappedPRB;
reservation.Period = ssburst.SSBPeriodicity * (gnb.SubcarrierSpacing/15); % Period in slots
pdsch.Reserved(end+1) = reservation;

%% Processing Loop
% To determine the throughput at each SNR point, the PDSCH data is analyzed 
% per transmission instance using the following steps:
% 
% * _Update current HARQ process._ Check the CRC of the previous
% transmission for the given HARQ process. Determine whether a
% retransmission is required. If that is not the case generate new data.
% * _Resource grid generation._ Channel coding is performed by
% <matlab:doc('nrDLSCH') nrDLSCH>. It operates on the input transport
% block provided. Internally, it keeps a copy of the transport block in
% case a retransmission is required. The coded bits are modulated on the 
% PDSCH by <matlab:doc('nrPDSCH') nrPDSCH>. The precoding operation is 
% applied to the resulting signal.
% * _Waveform generation._ The generated grid is then OFDM modulated.
% * _Noisy channel modeling._ The waveform is passed through a CDL or TDL
% fading channel. AWGN is added. The SNR is defined per RE at each UE
% antenna. For an SNR of 0 dB the signal and noise contribute equally to
% the energy per PDSCH RE per receive antenna [ <#14 5> ].
% * _Perform synchronization and OFDM demodulation._ For perfect
% synchronization, the channel impulse response is reconstructed and used
% to synchronize the received waveform. For practical synchronization, the
% received waveform is correlated with the PDSCH DM-RS. The synchronized
% signal is then OFDM demodulated.
% * _Perform channel estimation._ If perfect channel estimation is used,
% the channel impulse response is reconstructed and OFDM demodulated to
% provide a channel estimate. For practical channel estimation, the
% PDSCH DM-RS is used.
% * _Precoding matrix calculation._ The precoding matrix W for the next
% transmission is calculated using singular value decomposition (SVD). A
% single matrix is obtained for the full allocation by averaging the
% channel conditions. Therefore, for a channel with frequency selectivity,
% W could be less accurate for larger allocated bandwidths.
% * _Decode the PDSCH._ The recovered PDSCH symbols for all transmit and
% receive antenna pairs, along with a noise estimate, are demodulated and
% descrambled by <matlab:doc('nrPDSCHDecode') nrPDSCHDecode> to obtain an
% estimate of the received codewords.
% * _Decode the Downlink Shared Channel (DL-SCH) and store the block CRC
% error for a HARQ process._ The vector of decoded soft bits is passed to
% <matlab:doc('nrDLSCHDecoder') nrDLSCHDecoder> which decodes the codeword
% and returns the block CRC error used to determine the throughput of the
% system.

% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(snrIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(snrIn),1);

% Set up Redundancy Version (RV) sequence, number of HARQ processes and
% the sequence in which the HARQ processes are used
if pdsch.EnableHARQ
    % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
    % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
    % should be used. If self-decodability is the priority, it should be
    % taken into account that the upper limit of the code rate at which
    % each RV is self-decodable is in the following order: 0>3>>2>1
    rvSeq = [0 2 3 1];
else
    % HARQ disabled - single transmission with RV=0, no retransmissions
    rvSeq = 0; 
end

% Create DLSCH encoder system object
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = pdsch.TargetCodeRate;

% Create DLSCH decoder system object
% Use layered belief propagation for LDPC decoding, with half the number of
% iterations as compared to the default for belief propagation decoding
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = pdsch.TargetCodeRate;
%decodeDLSCH.LDPCDecodingAlgorithm = 'Layered belief propagation';
decodeDLSCH.LDPCDecodingAlgorithm = 'Normalized min-sum';
decodeDLSCH.MaximumLDPCIterationCount = 15;

% The temporary variables 'gnb_init', 'pdsch_init', 'ssburst_init' and
% 'decodeDLSCH_init' are used to create the temporary variables 'gnb',
% 'pdsch', 'ssburst' and 'decodeDLSCH' within the SNR loop to create
% independent instances in case of parallel simulation
gnb_init = gnb;
pdsch_init = pdsch;
ssburst_init = ssburst;
decodeDLSCH_init = clone(decodeDLSCH);

for snrIdx = 1:numel(snrIn) % comment out for parallel computing
% parfor snrIdx = 1:numel(snrIn) % uncomment for parallel computing
% To reduce the total simulation time, you can execute this loop in
% parallel by using the Parallel Computing Toolbox. Comment out the 'for'
% statement and uncomment the 'parfor' statement. If the Parallel Computing
% Toolbox is not installed, 'parfor' defaults to normal 'for' statement

    % Set the random number generator settings to default values
    rng('default');
    
    % Initialize variables for this SNR point, required for initialization
    % of variables when using the Parallel Computing Toolbox
    gnb = gnb_init;
    pdsch = pdsch_init;
    ssburst = ssburst_init;
    decodeDLSCH = clone(decodeDLSCH_init);
    pathFilters = [];
    ssbWaveform = [];
    
    SNRdB = snrIn(snrIdx);
    fprintf('\nSimulating transmission scheme 1 (%dx%d) and SCS=%dkHz with %s channel at %gdB SNR for %d 10ms frame(s)\n',...
        nTxAnts,nRxAnts,gnb.SubcarrierSpacing, ...
        channelType,SNRdB,gnb.NFrames); 
        
    % Initialize variables used in the simulation and analysis
    bitTput = [];           % Number of successfully received bits per transmission
    txedTrBlkSizes = [];    % Number of transmitted info bits per transmission
    
    % Specify the order in which we cycle through the HARQ processes
    NHARQProcesses = 16;
    harqSequence = 1:NHARQProcesses;

    % Initialize the state of all HARQ processes
    harqProcesses = hNewHARQProcesses(NHARQProcesses,rvSeq,gnb.NumCW);
    harqProcCntr = 0; % HARQ process counter
        
    % Reset the channel so that each SNR point will experience the same
    % channel realization
    reset(channel);
    
    % Total number of OFDM symbols in the simulation period
    waveformInfo = hOFDMInfo(gnb);
    NSymbols = gnb.NFrames * 10 * waveformInfo.SymbolsPerSubframe;
    
    % OFDM symbol number associated with start of each PDSCH transmission
    gnb.NSymbol = 0;
    
    % Running counter of the number of PDSCH transmission instances
    % The simulation will use this counter as the slot number for each
    % PDSCH
    pdsch.NSlot = 0;
    
    % Index to the start of the current set of SS burst samples to be
    % transmitted
    ssbSampleIndex = 1;
    
    % Obtain a precoding matrix (wtx) to be used in the transmission of the
    % first transport block
    estChannelGrid = getInitialChannelEstimate(gnb,nTxAnts,channel);    
    newWtx = getPrecodingMatrix(pdsch.PRBSet,pdsch.NLayers,estChannelGrid);
    
    % Timing offset, updated in every slot for perfect synchronization and
    % when the correlation is strong for practical synchronization
    offset = 0;
    
    while gnb.NSymbol < NSymbols  % Move to next slot, gnb.NSymbol increased in steps of one slot               
        
        % Generate a new SS burst when necessary
        if (ssbSampleIndex==1)        
            nSubframe = gnb.NSymbol / waveformInfo.SymbolsPerSubframe;
            ssburst.NFrame = floor(nSubframe / 10);
            ssburst.NHalfFrame = mod(nSubframe / 5,2);
            [ssbWaveform,~,ssbInfo] = hSSBurst(ssburst);
        end
        
        % Get HARQ process index for the current PDSCH from HARQ index table
        harqProcIdx = harqSequence(mod(harqProcCntr,length(harqSequence))+1);
        
        % Update current HARQ process information (this updates the RV
        % depending on CRC pass or fail in the previous transmission for
        % this HARQ process)
        harqProcesses(harqProcIdx) = hUpdateHARQProcess(harqProcesses(harqProcIdx),gnb.NumCW);
        
        % Calculate the transport block sizes for the codewords in the slot
        [pdschIndices,dmrsIndices,dmrsSymbols,pdschIndicesInfo] = hPDSCHResources(gnb,pdsch);
        trBlkSizes = hPDSCHTBS(pdsch,pdschIndicesInfo.NREPerPRB-Xoh_PDSCH);
        
        % HARQ: check CRC from previous transmission per codeword, i.e. is
        % a retransmission required?
        for cwIdx = 1:gnb.NumCW
            NDI = false;
            if harqProcesses(harqProcIdx).blkerr(cwIdx) % Errored
                if (harqProcesses(harqProcIdx).RVIdx(cwIdx)==1) % end of rvSeq
                    resetSoftBuffer(decodeDLSCH,cwIdx-1,harqProcIdx-1);
                    NDI = true;
                end
            else    % No error
                NDI = true;
            end
            if NDI 
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqProcIdx-1);
            end
        end
                                
        % Encode the DL-SCH transport blocks
        codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NLayers,...
            pdschIndicesInfo.G,harqProcesses(harqProcIdx).RV,harqProcIdx-1);

        % Get wtx (precoding matrix) calculated in previous slot
        wtx = newWtx;
        
        % PDSCH modulation and precoding
        pdschSymbols = nrPDSCH(codedTrBlock,pdsch.Modulation,pdsch.NLayers,gnb.NCellID,pdsch.RNTI);
        pdschSymbols = pdschSymbols*wtx;
        
        % PDSCH mapping in grid associated with PDSCH transmission period
        pdschGrid = zeros(waveformInfo.NSubcarriers,waveformInfo.SymbolsPerSlot,nTxAnts);
        [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
        pdschGrid(pdschAntIndices) = pdschSymbols;
        
        % PDSCH DM-RS precoding and mapping
        for p = 1:size(dmrsSymbols,2)
            [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
            pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols(:,p)*wtx(p,:);
        end
      
        % OFDM modulation of associated resource elements
        txWaveform = hOFDMModulate(gnb, pdschGrid);
        
        % Add the appropriate portion of SS burst waveform to the
        % transmitted waveform
        Nt = size(txWaveform,1);
        txWaveform = txWaveform + ssbWaveform(ssbSampleIndex + (0:Nt-1),:);
        ssbSampleIndex = mod(ssbSampleIndex + Nt,size(ssbWaveform,1));

        % Pass data through channel model. Append zeros at the end of the
        % transmitted waveform to flush channel content. These zeros take
        % into account any delay introduced in the channel. This is a mix
        % of multipath delay and implementation delay. This value may 
        % change depending on the sampling rate, delay profile and delay
        % spread
        txWaveform = [txWaveform; zeros(maxChDelay, size(txWaveform,2))];
        [rxWaveform,pathGains,sampleTimes] = channel(txWaveform);
        
        % Add AWGN to the received time domain waveform
        % Normalize noise power to take account of sampling rate, which is
        % a function of the IFFT size used in OFDM modulation. The SNR
        % is defined per RE for each receive antenna (TS 38.101-4).
        SNR = 10^(SNRdB/20);    % Calculate linear noise gain
        N0 = 1/(sqrt(2.0*nRxAnts*double(waveformInfo.Nfft))*SNR);
        noise = N0*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
        rxWaveform = rxWaveform + noise;
       
        if (perfectChannelEstimator)
            % Perfect synchronization. Use information provided by the channel
            % to find the strongest multipath component
            pathFilters = getPathFilters(channel); % get path filters for perfect channel estimation
            [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
            
        else
            % Practical synchronization. Correlate the received waveform 
            % with the PDSCH DM-RS to give timing offset estimate 't' and 
            % correlation magnitude 'mag'. The function
            % hSkipWeakTimingOffset is used to update the receiver timing
            % offset. If the correlation peak in 'mag' is weak, the current
            % timing estimate 't' is ignored and the previous estimate
            % 'offset' is used
            [t,mag] = nrTimingEstimate(rxWaveform,gnb.NRB,gnb.SubcarrierSpacing,pdsch.NSlot,dmrsIndices,dmrsSymbols,'CyclicPrefix',gnb.CyclicPrefix); 
            offset = hSkipWeakTimingOffset(offset,t,mag);
                    
        end

        
%        collect_data=1;
%        global H_idl;
%        global H_hard_wnd;
       if(collect_data)
           pathFilters = getPathFilters(channel); % get path filters for perfect channel estimation
            [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
            rxWaveform = rxWaveform(1+offset:end, :);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid, including padding in the event that practical
        % synchronization results in an incomplete slot being demodulated
        rxGrid = hOFDMDemodulate(gnb, rxWaveform);
        [K,L,R] = size(rxGrid);
        if (L < waveformInfo.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,waveformInfo.SymbolsPerSlot-L,R));
        end
       H_idl_tmp = nrPerfectChannelEstimate(pathGains,pathFilters,gnb.NRB,gnb.SubcarrierSpacing,pdsch.NSlot,offset,sampleTimes,gnb.CyclicPrefix);
      
       [t,mag] = nrTimingEstimate(rxWaveform,gnb.NRB,gnb.SubcarrierSpacing,pdsch.NSlot,dmrsIndices,dmrsSymbols,'CyclicPrefix',gnb.CyclicPrefix); %#ok<UNRCH>
            offset = hSkipWeakTimingOffset(offset,t,mag);
                    rxWaveform = rxWaveform(1+offset:end, :);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid, including padding in the event that practical
        % synchronization results in an incomplete slot being demodulated
        rxGrid = hOFDMDemodulate(gnb, rxWaveform);
        [K,L,R] = size(rxGrid);
        if (L < waveformInfo.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,waveformInfo.SymbolsPerSlot-L,R));
        end
       [H_hard_wnd_tmp,~] = nrChannelEstimate(rxGrid,dmrsIndices,dmrsSymbols,'CyclicPrefix',gnb.CyclicPrefix,'CDMLengths',pdschIndicesInfo.CDMLengths); %#ok<UNRCH>
       
       %visualization of channel estimation from first antenna
%        if(1)
%            
%        figure;
%        subplot(1,2,1);
%        imagesc(abs(H_hard_wnd_tmp(:,:,1)));
%        xlabel('OFDM Symbol');
%        ylabel('Subcarrier');
%        title('Hard Window Estimate Magnitude');
%        subplot(1,2,2);
%        imagesc(abs(H_idl_tmp(:,:,1)));
%        xlabel('TTIs');
%        ylabel('Subcarrier');
%        title('Ideal Estimate Magnkitude');
%        end
       H_idl=[H_idl H_idl_tmp];
       H_hard_wnd=[H_hard_wnd H_hard_wnd_tmp];
        %data for 2D NN, subcariers by antenna 
     
%               if(1)
%            
%        figure;
%        subplot(1,2,1);
%        imagesc(abs(H_hard_wnd2d(:,:)));
%        xlabel('OFDM Symbol');
%        ylabel('Subcarrier');
%        title('Hard Window Estimate Magnitude');
%        subplot(1,2,2);
%        imagesc(abs(H_idl2d(:,:)));
%        xlabel('TTIs');
%        ylabel('Subcarrier');
%        title('Ideal Estimate Magnkitude');
%        end
       
        
       
       
           
         if mod(size(H_idl,2)/14,save_period)==0
           save ('H_store_minus5dB_doppler30_'+ string(size(H_idl,2)/14)+  '.mat', 'H_idl','H_hard_wnd');
            if(collect_data2d)
             H_hard_wnd2d(1:size(H_hard_wnd,1),1:size(H_hard_wnd,3)*size(H_hard_wnd,2))=zeros();
             H_idl2d(1:size(H_idl,1),1:size(H_idl,3)*size(H_idl,2))=zeros();
              for i=1:size(H_idl,2)
                H_idl2d(1:size(H_idl,1),4*(i-1)+1:4*i)=squeeze(H_idl(:,i,:));
                H_hard_wnd2d(1:size(H_idl,1),4*(i-1)+1:4*i)=squeeze(H_hard_wnd(:,i,:));
              end
              save ('2DH_store_minus5dB_doppler30'+ string(size(H_idl,2)/14)+  '.mat', 'H_idl2d','H_hard_wnd2d');
            end
            
         end
       else
       rxWaveform = rxWaveform(1+offset:end, :);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid, including padding in the event that practical
        % synchronization results in an incomplete slot being demodulated
        rxGrid = hOFDMDemodulate(gnb, rxWaveform);
        [K,L,R] = size(rxGrid);
           if (L < waveformInfo.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,waveformInfo.SymbolsPerSlot-L,R));
           end    
       end
        if   (perfectChannelEstimator)
            % Perfect channel estimation, using the value of the path gains
            % provided by the channel. This channel estimate does not
            % include the effect of transmitter precoding
            estChannelGrid = nrPerfectChannelEstimate(pathGains,pathFilters,gnb.NRB,gnb.SubcarrierSpacing,pdsch.NSlot,offset,sampleTimes,gnb.CyclicPrefix);

            % Get perfect noise estimate (from the noise realization)
            noiseGrid = hOFDMDemodulate(gnb,noise(1+offset:end ,:));
            noiseEst = var(noiseGrid(:));

            % Get precoding matrix for next slot
            newWtx = getPrecodingMatrix(pdsch.PRBSet,pdsch.NLayers,estChannelGrid);

            % Apply precoding to estChannelGrid
            estChannelGrid = precodeChannelEstimate(estChannelGrid,wtx.');
        else
            % Practical channel estimation between the received grid and
            % each transmission layer, using the PDSCH DM-RS for each
            % layer. This channel estimate includes the effect of
            % transmitter precoding
            if (perfect_noise_estimation) 
            % Use perfect noise estimate (from the noise realization)
            noiseGrid = hOFDMDemodulate(gnb,noise(1+offset:end ,:));
            noiseEst = var(noiseGrid(:));
            else
            %Use practical noise estimation
            [estChannelGrid,noiseEst] = nrChannelEstimate(rxGrid,dmrsIndices,dmrsSymbols,'CyclicPrefix',gnb.CyclicPrefix,'CDMLengths',pdschIndicesInfo.CDMLengths); %#ok<UNRCH>
           
            end
            
            %Array of size Subcariers (RB*12) x TTIs (OFDM symbols 14 TTis)
            %x Rx antennas
          
            if (inference)


            loadtrainedDnCNN.mat;
            RB=20;
            RX=4;
            for i=1:size(estChannelGrid,2)/14
            	H_hard_wnd_train_real(:,1:28,1:RX,1,i)=[real( estChannelGrid(:,:,1:RX,i)) imag( estChannelGrid(:,:,1:RX,i))] ; 
            end   
            
            Ypred = predict(net,H_hard_wnd_train_real);    
            H_NN=estChannelGrid-double(Ypred(1:RB*12,1:14,1:RX))-double(1i*Ypred(1:RB*12,15:28,1:RX));
            %estChannelGrid=H_NN;
            estChannelGrid=H_NN;
       %visualization of channel estimations from first antenna
%        if(0)
%            
%        figure;
%        subplot(2,2,1);
%        imagesc(abs(estChannelGrid(:,:,1)));
%        xlabel('OFDM Symbol');
%        ylabel('Subcarrier');
%        title('Hard Window Estimate Magnitude');
%        subplot(2,2,2);
%        imagesc(abs(H_NN(:,:,1)));
%        xlabel('TTIs');
%        ylabel('Subcarrier');
%        title('NN Estimate Magnitude');
%        subplot(2,2,[3,4]);
%        
%        imagesc(abs(H_idl_tmp(:,:,1)));
%        xlabel('TTIs');
%        ylabel('Subcarrier');
%        title('Ideal Estimate Magnitude');
%        end
            
            elseif  (inference2d)

               load trainedDnCNN.mat;
               
               
               RB=20;
               RX=4;
               
%               for i=1:size(estChannelGrid,2)
%                 H_hard_wnd2d(1:size(estChannelGrid,1),4*(i-1)+1:4*i)=squeeze(estChannelGrid(:,i,:));
%               end
% 
%                 H_hard_wnd_train_real(1:RB*12,1:RX)=real(H_hard_wnd2d);
%                 H_hard_wnd_train_real(1:RB*12,RX+1:RX*2)=imag(H_hard_wnd2d);
                H_NN_grid(1:RB*12,1:size(estChannelGrid,2),1:RX)=zeros();
                   for i=1:size(estChannelGrid,2)
                H_hard_wnd2dinference(1:size(estChannelGrid,1),1:size(estChannelGrid,3))=squeeze(estChannelGrid(:,i,:));
                H_hard_wnd2dinference_real(1:RB*12,1:RX)=real(H_hard_wnd2dinference);
                H_hard_wnd2dinference_real(1:RB*12,RX+1:RX*2)=imag(H_hard_wnd2dinference);
                Ypred = predict(net,H_hard_wnd2dinference_real);  
                H_NN_grid(1:RB*12,i,1:RX)=H_hard_wnd2dinference-double(Ypred(1:RB*12,1:RX,1,1))-double(1i*Ypred(1:RB*12,RX+1:RX*2,1,1));
     
             
               
                   end

                 estChannelGrid=H_NN_grid;

                 
             
            
            
                
            
        else
            [estChannelGrid,noiseEst] = nrChannelEstimate(rxGrid,dmrsIndices,dmrsSymbols,'CyclicPrefix',gnb.CyclicPrefix,'CDMLengths',pdschIndicesInfo.CDMLengths); %#ok<UNRCH>
            
            
            
            
            
            
            % Remove precoding from estChannelGrid prior to precoding
            % matrix calculation
            estChannelGridPorts = precodeChannelEstimate(estChannelGrid,conj(wtx));
            
            % Get precoding matrix for next slot
            newWtx = getPrecodingMatrix(pdsch.PRBSet,pdsch.NLayers,estChannelGridPorts);
        end
       end
        % Get PDSCH resource elements from the received grid
        [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChannelGrid);
        
        % Equalization
        [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
        
        % Decode PDSCH physical channel
        [dlschLLRs,rxSymbols] = nrPDSCHDecode(pdschEq,pdsch.Modulation,gnb.NCellID,pdsch.RNTI,noiseEst);

        % Scale LLRs by CSI
        csi = nrLayerDemap(csi); % CSI layer demapping
        for cwIdx = 1:gnb.NumCW
            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % bits per symbol
            csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);   % expand by each bit per symbol
            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % scale
        end
        
        % Decode the DL-SCH transport channel
        decodeDLSCH.TransportBlockLength = trBlkSizes;
        [decbits,harqProcesses(harqProcIdx).blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NLayers,harqProcesses(harqProcIdx).RV,harqProcIdx-1);
        
        % Store values to calculate throughput (only for active PDSCH instances)
        if(any(trBlkSizes ~= 0))
            bitTput = [bitTput trBlkSizes.*(1-harqProcesses(harqProcIdx).blkerr)];
            txedTrBlkSizes = [txedTrBlkSizes trBlkSizes];                         
        end
        
        % Update starting symbol number of next PDSCH transmission
        gnb.NSymbol = gnb.NSymbol + size(pdschGrid,2);
        % Update count of overall number of PDSCH transmissions
        pdsch.NSlot = pdsch.NSlot + 1;
        % Update HARQ process counter
        harqProcCntr = harqProcCntr + 1;
                
        % Display transport block error information per codeword managed by current HARQ process
        fprintf('\n(%3.2f%%) HARQ Proc %d: ',100*gnb.NSymbol/NSymbols,harqProcIdx);
        estrings = {'passed','failed'};
        rvi = harqProcesses(harqProcIdx).RVIdx; 
        for cw=1:length(rvi)
            cwrvi = rvi(cw);
            % Create a report on the RV state given position in RV sequence and decoding error
            if cwrvi == 1
                ts = sprintf('Initial transmission (RV=%d)',rvSeq(cwrvi));
            else
                ts = sprintf('Retransmission #%d (RV=%d)',cwrvi-1,rvSeq(cwrvi));
            end
            fprintf('CW%d:%s %s. ',cw-1,ts,estrings{1+harqProcesses(harqProcIdx).blkerr(cw)}); 
        end
         
     end
    
    % Calculate maximum and simulated throughput
    maxThroughput(snrIdx) = sum(txedTrBlkSizes); % Max possible throughput
    simThroughput(snrIdx) = sum(bitTput,2);      % Simulated throughput
    
    % Display the results dynamically in the command window
    fprintf([['\n\nThroughput(Mbps) for ', num2str(gnb.NFrames) ' frame(s) '],...
        '= %.4f\n'], 1e-6*simThroughput(snrIdx)/(gnb.NFrames*10e-3));
    fprintf(['Throughput(%%) for ', num2str(gnb.NFrames) ' frame(s) = %.4f\n'],...
        simThroughput(snrIdx)*100/maxThroughput(snrIdx));
     
end

%% Results
% Display the measured throughput. This is calculated as the percentage of
% the maximum possible throughput of the link given the available resources
% for data transmission.

figure;
plot(snrIn,simThroughput*100./maxThroughput,'o-.')
xlabel('SNR (dB)'); ylabel('Throughput (%)'); grid on;
title(sprintf('(%dx%d) / NRB=%d / SCS=%dkHz',...
              nTxAnts,nRxAnts,gnb_init.NRB,gnb_init.SubcarrierSpacing));

% Bundle key parameters and results into a combined structure for recording
simResults.simParameters = simParameters;
simResults.simThroughput = simThroughput;

%%
% The figure below shows throughput results obtained simulating 10000
% subframes (|NFrames = 1000|, |SNRIn = -18:2:16|).
%
% <<../longRunThroughput.png>>
%

%% Appendix
% This example uses the following helper functions:
% 
% * <matlab:edit('hArrayGeometry.m') hArrayGeometry.m>
% * <matlab:edit('hNewHARQProcesses.m') hNewHARQProcesses.m>
% * <matlab:edit('hOFDMDemodulate.m') hOFDMDemodulate.m>
% * <matlab:edit('hOFDMInfo.m') hOFDMInfo.m>
% * <matlab:edit('hOFDMModulate.m') hOFDMModulate.m>
% * <matlab:edit('hPDSCHResources.m') hPDSCHResources.m>
% * <matlab:edit('hPDSCHTBS.m') hPDSCHTBS.m>
% * <matlab:edit('hSkipWeakTimingOffset.m') hSkipWeakTimingOffset.m>
% * <matlab:edit('hSSBurst.m') hSSBurst.m>
% * <matlab:edit('hUpdateHARQProcess.m') hUpdateHARQProcess.m>

%% Selected Bibliography
% # 3GPP TS 38.211. "NR; Physical channels and modulation (Release 15)."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.212. "NR; Multiplexing and channel coding (Release 15)." 3rd
% Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.213. "NR; Physical layer procedures for control (Release
% 15)." 3rd Generation Partnership Project; Technical Specification Group
% Radio Access Network.
% # 3GPP TS 38.214. "NR; Physical layer procedures for data (Release 15)."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # R1-166999. "Detailed configuration of F-OFDM and W-OFDM for LLS
% evaluation", 3GPP RAN WG1 #86, Spreadtrum Communications, August 2016.
% # 3GPP TR 38.901. "Study on channel model for frequencies from 0.5 to 100
% GHz (Release 15)." 3rd Generation Partnership Project; Technical
% Specification Group Radio Access Network.
% # 3GPP TS 38.101-4. "NR; User Equipment (UE) radio transmission and
% reception. Part 4: Performance requirements (Release 15)." 3rd Generation
% Partnership Project; Technical Specification Group Radio Access Network.

%% Local Functions

function validateNLayers(simParameters)
% Validate the number of layers
    if length(simParameters.PDSCH.PortSet)~= simParameters.PDSCH.NLayers
        error('The number of elements of PortSet (%d) must be the same as the number of layers (%d)',...
            length(simParameters.PDSCH.PortSet), simParameters.PDSCH.NLayers);
    end

    if simParameters.PDSCH.NLayers > min(simParameters.NTxAnts,simParameters.NRxAnts)
        error('The number of layers (%d) must satisfy NLayers <= min(NTxAnts,NRxAnts) = min(%d,%d) = (%d)',...
            simParameters.PDSCH.NLayers,simParameters.NTxAnts,simParameters.NRxAnts,min(simParameters.NTxAnts,simParameters.NRxAnts));
    end
end

function estChannelGrid = getInitialChannelEstimate(gnb,nTxAnts,channel)
% Obtain channel estimate before first transmission. This can be used to
% obtain a precoding matrix for the first slot.

    ofdmInfo = hOFDMInfo(gnb);
    
    chInfo = info(channel);
    maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
    
    % Temporary waveform (only needed for the sizes)
    tmpWaveform = zeros((ofdmInfo.SamplesPerSubframe/ofdmInfo.SlotsPerSubframe)+maxChDelay,nTxAnts);
    
    % Filter through channel    
    [~,pathGains,sampleTimes] = channel(tmpWaveform);
    
    % Perfect timing synch    
    pathFilters = getPathFilters(channel);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    nslot = gnb.NSymbol/ofdmInfo.SymbolsPerSlot;
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(pathGains,pathFilters,gnb.NRB,gnb.SubcarrierSpacing,nslot,offset,sampleTimes);
    
end

function wtx = getPrecodingMatrix(PRBSet,NLayers,hestGrid)
% Calculate precoding matrix given an allocation and a channel estimate
    
    % Allocated subcarrier indices
    allocSc = (1:12)' + 12*PRBSet(:).';
    allocSc = allocSc(:);
    
    % Average channel estimate
    [~,~,R,P] = size(hestGrid);
    estAllocGrid = hestGrid(allocSc,:,:,:);
    Hest = permute(mean(reshape(estAllocGrid,[],R,P)),[2 3 1]);
    
    % SVD decomposition
    [~,~,V] = svd(Hest);
    wtx = V(:,1:NLayers).';
    wtx = wtx / sqrt(NLayers); % Normalise by NLayers

end

function estChannelGrid = precodeChannelEstimate(estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate

    % Linearize 4D matrix and reshape after multiplication
    K = size(estChannelGrid,1);
    L = size(estChannelGrid,2);
    R = size(estChannelGrid,3);
    estChannelGrid = reshape(estChannelGrid,K*L*R,[]);
    estChannelGrid = estChannelGrid * W;
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end

function [mappedPRB,mappedSymbols] = mapNumerology(subcarriers,symbols,nrbs,nrbt,fs,ft)
% Map the SSBurst numerology to PDSCH numerology. The outputs are:
%   - mappedPRB: 0-based PRB indices for carrier resource grid (arranged in a column)
%   - mappedSymbols: 0-based OFDM symbol indices in a slot for carrier resource grid (arranged in a row)
%     carrier resource grid is sized using gnb.NRB, gnb.CyclicPrefix, spanning 1 slot
% The input parameters are:
%   - subcarriers: 1-based row subscripts for SSB resource grid (arranged in a column)
%   - symbols: 1-based column subscripts for SSB resource grid (arranged in an N-by-4 matrix, 4 symbols for each transmitted burst in a row, N transmitted bursts)
%     SSB resource grid is sized using ssbInfo.NRB, normal CP, spanning 5 subframes
%   - nrbs: source (SSB) NRB
%   - nrbt: target (carrier) NRB
%   - fs: source (SSB) SCS
%   - ft: target (carrier) SCS

    mappedPRB = unique(fix((subcarriers-(nrbs*6) - 1)*fs/(ft*12) + nrbt/2),'stable');
    
    symbols = symbols.';
    symbols = symbols(:).' - 1;

    if (ft < fs)
        % If ft/fs < 1, reduction
        mappedSymbols = unique(fix(symbols*ft/fs),'stable');
    else
        % Else, repetition by ft/fs
        mappedSymbols = reshape((0:(ft/fs-1))' + symbols(:)'*ft/fs,1,[]);
    end
    
end

