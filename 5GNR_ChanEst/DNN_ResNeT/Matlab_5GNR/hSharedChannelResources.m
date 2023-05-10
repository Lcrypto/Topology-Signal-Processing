%hSharedChannelResources 5G NR PDSCH/PUSCH, DM-RS and PT-RS resource element indices, DM-RS and PT-RS values
%   [IND,DMRSIND,DMRS,PTRSIND,PTRS,INFO,GINFO] = hSharedChannelResources(CFT,GNB,CHS)
%   returns the resource element (RE) indices for a 5G NR PDSCH/PUSCH
%   transmission, associated DM-RS and PT-RS, given the time (symbols) and
%   frequency (PRBs) allocation of the channels, the DM-RS configuration
%   and the PT-RS configuration. This function implements aspects common to
%   both PDSCH and PUSCH, and supports OFDM based transmission without any
%   uplink extensions to cover SC-FDMA or frequency hopping. It is used by
%   hPDSCHResources and hPUSCHResources to provide the overall PDSCH and
%   PUSCH resourcing, DM-RS and PT-RS functionality.
%
%   The 1-based linear indices are returned in matrix IND. They are
%   defined relative to a three-dimensional RE grid representing a 
%   14/12-symbol slot for the full carrier or bandwidth (in the PDSCH/PUSCH
%   numerology) across the layers/DM-RS ports of the transmission. Each
%   column of IND represents the grid locations for a separate layer/port
%   (the third dimension of the grid). The DM-RS RE indices have the same
%   format and are returned in matrix DMRSIND. The complex values of DM-RS
%   sequence are also returned in matrix DMRS. 
%
%   The callback configuration input, CFT, must be a structure including
%   the fields:
%   MappingType      - Slot mapping parameter name 
%                      ('PDSCHMappingType','PUSCHMappingType')
%   dmrsSymbolsTable - Function handle to provide l^bar DM-RS position indices
%   getDMRSSequence  - Optional. Function handle to override the generation
%                      of the full set of DM-RS sequence values
%                      (if parameter field is not provided or if empty then
%                      PDSCH/PUSCH PRBS based DM-RS for CP-OFDM are created)
%   getPTRSResources - Optional. Function handle to override the generation
%                      of PT-RS symbols and indices
%                      (if parameter field is not provided or if empty then
%                      PDSCH/PUSCH PT-RS symbols and indices for CP-OFDM
%                      are created, provided EnablePTRS is set to 1)
%
%   The general settings input, GNB, must be a structure including
%   the fields:
%   NRB             -   Number of resource blocks for the carrier
%                       or bandwidth part (in the PDSCH/PUSCH numerology)
%   CyclicPrefix      - Optional. Cyclic prefix length 
%                       ('Normal'(default),'Extended')
%   SubcarrierSpacing - Optional. Subcarrier Spacing (kHz)
%                       (15(default),30,60,120,240)
%   RBOffset          - Optional. Position of BWP in SCS carrier (default 0)
%                       (Only used if CHS.VRBToPRBInterleaving is defined
%                       and true)
%
%   The channel specific input, CHS, must be a structure including the
%   fields:
%   NSlot                   - Optional. Slot number of transmission (default 0) 
%   PRBSet                  - PRBs allocated to the channel (0-based indices)
%   PRBRefPoint             - Optional. PRB index that PRBSet is relative to
%                             (0-based) (default 0)
%   SymbolSet               - OFDM symbols allocated to transmission within slot
%                             (including DM-RS symbols, 0-based indices, max range 0...13)
%   PortSet                 - DM-RS ports used by PDSCH/PUSCH
%                             (0-based indices, max range 0...11,
%                             mapping to ports p=1000...1011/p=0...11 respectively)
%   NLayers                 - Optional. Number of layers
%                             (Only used if PortSet has not been defined)
%   Modulation              - Modulation type(s) of codeword(s)
%                             ('pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM')
%   Reserved                - Optional. Reserved PRB patterns
%                             (structure array, see below)
%   DMRSConfigurationType   - DM-RS configuration type (1,2)
%   NumCDMGroupsWithoutData - Optional. Number of CDM groups without data
%                             (0(default),1,2,3)
%   VRBToPRBInterleaving    - Optional. Enable or disable interleaved
%                             resource mapping (false(default), true)
%   VRBBundleSize           - Optional. vrb-ToPRB-Interleaver parameter (2(default), 4)
%                             (Only used if CHS.VRBToPRBInterleaving is true)
%   RNTI                    - Optional. Radio network temporary identifier
%                             (0...65535) (default 0)
%   TransformPrecoding      - Optional. Transform precoding flag (0(default),1)
%   IntraSlotFreqHopping    - Optional. Intra-slot frequency hopping
%                             ('disabled'(default),'enabled')
%
%   If CFT.getDMRSSequence does not override the generation of the DM-RS
%   values then the following scrambling identity parameters are required:
%   NIDNSCID                - DM-RS scrambling identity (0...65535)
%   NSCID                   - DM-RS scrambling initialization (0,1)
%
%   The DM-RS OFDM symbol locations can either be defined explicitly 
%   using the single parameter:
%   DMRSSymbolSet           - OFDM symbols containing the DM-RS within the 
%                             shared channel allocation in slot 
%                             (0-based indices, max range 0...13)
%   or, defined implicitly via the FT.dmrsSymbolsTable function handle and 
%   its own parameter requirements, and the following DM-RS parameters: 
%   DMRSTypeAPosition       - Mapping type A only. First DM-RS symbol position
%                             (2,3)
%   DMRSLength              - Number of front-loaded DM-RS symbols
%                             (1(single symbol),2(double symbol))
%
%   The following fields specific to PT-RS must be provided in CHS input:
%   EnablePTRS              - Optional. Enable PT-RS (0(default),1)
%   PTRSTimeDensity         - Optional. Time density (L_PT-RS) of PT-RS (1(default),2,4)
%   PTRSFrequencyDensity    - Optional. Frequency density (K_PT-RS) of PT-RS (2(default),4)
%   PTRSNumSamples          - Optional. Number of PT-RS samples in a group
%                             (NGroupSamp) (2(default),4)
%                             (Used only when transform precoding is enabled)
%   PTRSNumGroups           - Optional. Number of PT-RS groups (NPTRSGroup) (2(default),4,8) 
%                             (Used only when transform precoding is enabled)
%   PTRSREOffset            - Optional. Resource element offset ('00'(default),'01','10','11')
%   PTRSPortSet             - Optional. PT-RS antenna ports. It must be a
%                             subset of DM-RS antenna ports
%                             (default is lowest DM-RS port number)
%   PTRSNID                 - Optional. PT-RS scrambling identity (0...1007)
%                             (default is the DM-RS scrambling identity NRSID)
%                             (Used only when transform precoding is enabled)
%
%   Periodically recurring patterns of reserved PRB can defined using the
%   'Reserved' parameter. These PRB will be excluded from the generated 
%   indices and the transport and physical processing should rate-match 
%   around them. This parameter takes the format of a structure array where
%   each element defines a separate pattern. Each element, i, of the 
%   array should contain the following fields:
%   Reserved(i).PRB     - Reserved PRB (0-based indices, defined as a 
%                         vector or cell array)
%   Reserved(i).Symbols - OFDM symbols associated with reserved PRB 
%                         (0-based indices, spanning one or more slots)
%   Reserved(i).Period  - Total number of slots in the pattern period
% 
%   The reserved PRB indices can be specified as a vector or a cell array.
%   If a vector then the same PRBs are excluded in each OFDM symbol in the 
%   pattern. If the PRB indices are defined as a cell array then each cell
%   specifies the excluded PRBs for the associated OFDM symbol in the 
%   pattern. In the latter case, the length of the PRB cell array should
%   match the length of the 'Symbols' field, i.e. an individual set of PRBs
%   is defined for each reserved symbol. The symbols that form the 
%   time-domain locations of the reserved pattern can be greater than 13
%   and therefore cover multiple slots. The overall pattern will repeat
%   itself every 'Period' slots. If this field is empty then the pattern 
%   will not cyclically repeat itself.
% 
%   In terms of frequency domain DM-RS density, there are two different RRC
%   signaled configuration types ('dmrs_Type'). Configuration
%   type 1 defines 6 subcarriers per PRB per antenna port, comprising
%   alternating subcarriers. Configuration type 2 defines 4 subcarriers per
%   PRB per antenna ports, consisting of 2 groups of 2 neighboring
%   subcarriers. Different shifts are applied to the sets of subcarriers
%   used, dependent on the associated antenna port or CDM group. For
%   type 1, there are 2 possible CDM groups/shifts across up to 8 possible
%   antenna ports (p=0...11 for uplink/p=1000...1007 for downlink), and, 
%   for type 2, there are 3 possible CDM groups/shifts across 12 ports
%   (p=0...11/p=1000...1011). See TS 38.211 sections 6.4.1.1 and 7.4.1.1
%   for the full configuration details.
% 
%   INFO is the output structure containing the fields:
%   G             - Bit capacity of the PDSCH/PUSCH. This is should be the
%                   length of codeword to be output from the DL-SCH/UL-SCH 
%                   transport channel
%   Gd            - Number of resource elements per layer/port, equal to 
%                   the number of rows in the PDSCH/PUSCH indices
%   DMRSSymbolSet - The symbol numbers in a slot containing DM-RS (0-based)
%   NREPerPRB     - Number of RE per PRB allocated to PDSCH/PUSCH (not 
%                   accounting for any reserved resources)
%   CDMGroups     - CDM groups associated with the DM-RS antenna ports
%   CDMLengths    - A 2-element row vector [FD TD] specifying the length of
%                   FD-CDM and TD-CDM despreading required during channel
%                   estimation. The values depend on the frequency and
%                   time masks applied to groups of antenna ports
%   PTRSSymbolSet - The symbol numbers in a slot containing PT-RS (0-based)
%
%   GINFO is the output structure containing the fields:
%   SymbolsPerSlot - Number of OFDM symbols in a slot, 14 or 12 for normal
%                    and extended cyclic prefix, respectively
%
%   Example:
%   % Configure the function to create PRBS based DM-RS for 3 OFDM symbols 
%   % in a slot. Display the associated data channel and DM-RS resource
%   % element locations.
% 
%   % Define a callback configuration structure
%   ft = struct();
%   % Mapping type parameter name
%   ft.MappingType = 'PDSCHMappingType';    
%   % As a simple example, set up a function handle which will always
%   % returns 3 DM-RS position indices 0, 5 and 10.
%   ft.dmrsSymbolsTable = @(pxsch,typeb,nsymbols)[0 5 10];
%   % Don't override the default PRBS based DM-RS sequence generation
%   ft.getDMRSSequence = [];
%   % Don't override the PT-RS resources
%   ft.getPTRSResources = [];
%
%   % Set the number of resource blocks in carrier or BWP and 
%   % the numerology
%   gnb = struct('NRB',50);
%   gnb.SubcarrierSpacing = 15;
%   gnb.CyclicPrefix = 'Normal';
% 
%   % Get the number of OFDM symbols in a slot
%   symbperslot = sum(strcmpi(gnb.CyclicPrefix,["Normal","Extended"]) .* [14 12]);
%
%   % Specify the basic allocation properties to be full band, full
%   % slot using 2 layers/ports
%   pxsch = struct();
%   pxsch.NSlot = 0;                    % Slot number
%   pxsch.PRBSet = 0:gnb.NRB-1;         % Full band PRB allocation
%   pxsch.SymbolSet = 0:symbperslot-1;  % Full slot allocation
%   pxsch.PortSet = [0,2];              % Use DM-RS ports 1000 and 1002
%   pxsch.Modulation = 'QPSK';          % Modulation type
%
%   % Exclude some PRBs (e.g. for a CORESET in the downlink) in the middle
%   % of the carrier in the first two symbols in the slot
%   pxsch.Reserved.PRB = fix(gnb.NRB/2) + (-4:5);
%   pxsch.Reserved.Symbols = [0,1];
%   pxsch.Reserved.Period = 1;
% 
%   % Configure the DM-RS for config type 1 and slot-wise, type A
%   % mapping. Use double symbols for front-loaded DM-RS and configure an
%   % additional symbol pair towards the end of the slot
%   pxsch.DMRSConfigurationType = 1; % DM-RS configuration type (1,2)
%   pxsch.NIDNSCID = 1;            % DM-RS scrambling identity (0...65535)
%   pxsch.NSCID = 0;               % DM-RS scrambling initialization (0,1)
%   pxsch.PDSCHMappingType = 'A';  % Slot-wise mapping type (type A)
%   pxsch.DMRSTypeAPosition = 2;   % First DM-RS symbol position for type A 
%   pxsch.DMRSLength = 1;          % Specify single front-loaded DM-RS
%
%   % Configure the PT-RS with frequency density set to 2, time density set
%   % to 1, resource element offset set to '00' and PT-RS antenna port set
%   % to 0
%   pxsch.EnablePTRS = 1;           % Enable or disable PT-RS (1 or 0)  
%   pxsch.PTRSFrequencyDensity = 2; % Frequency density (2,4)
%   pxsch.PTRSTimeDensity = 1;      % Time density (1,2,4)
%   pxsch.PTRSREOffset = '00';      % Resource element offset('00','01','10','11')
%   pxsch.PTRSPortSet = 0;          % Antenna port set for PT-RS
%
%   % Display data, DM-RS and PT-RS RE locations of the first port of the
%   % grid
%   slotgrid = zeros(12*gnb.NRB,symbperslot,length(pxsch.PortSet));
%   [ind,dmrsind,dmrs,ptrsind,ptrs,info] = hSharedChannelResources(ft,gnb,pxsch);
%   slotgrid(ind) = 20;                 % Use light blue for data RE 
%   slotgrid(dmrsind) = 40*abs(dmrs);   % Use green for DM-RS RE
%   slotgrid(ptrsind) = 70*abs(ptrs);   % Use yellow for PT-RS RE
%   figure;
%   imagesc(slotgrid(:,:,1));
%   title('Data, DM-RS and PT-RS resource elements');
%   axis('xy'), xlabel('Symbols'), ylabel('Subcarriers');
%  
%   See also hPDSCHResources, hPUSCHResources.

%   Copyright 2019 The MathWorks, Inc.
 
function [pxschIndices,dmrsIndices,dmrSymbols,ptrsIndices,ptrsSymbols,pxschInfo,genInfo] ...
                                      = hSharedChannelResources(ftable,gnb,pxsch)

    % Capture the SCS and number of slots in a 10ms frame
    if isfield(gnb,'SubcarrierSpacing')
        scs = gnb.SubcarrierSpacing;
    else
        scs = 15;  % Default to 15kHz SCS
    end
    slotsperframe = 10*(scs/15);
    
    % Establish the number of OFDM symbols in a single slot
    if isfield(gnb,'CyclicPrefix')
        cpoptions = ["Normal","Extended"];
        cpselected = strcmpi(gnb.CyclicPrefix,cpoptions);
        if ~any(cpselected)
            error('The cyclic prefix (%s) must be one of (%s)' , string(gnb.CyclicPrefix), join(cpoptions,','));
        end
        symbperslot = sum(cpselected .* [14 12]);
    else
        symbperslot = 14;   % Default to normal CP
    end
    
    % Capture the set of active OFDM symbols (within the slot) associated
    % with the allocation
    if isfield(pxsch,'SymbolSet')
        usedsymbols = pxsch.SymbolSet < symbperslot;
        allocatedsymbols = pxsch.SymbolSet(usedsymbols);        
    else
        usedsymbols = ones(1,symbperslot);
        allocatedsymbols = (0:symbperslot-1);   % Default to a full slot allocation
        pxsch.SymbolSet = allocatedsymbols;
    end
   
    % Capture set of PRB (within the slot) associated with the allocation
    if iscell(pxsch.PRBSet)
        if length(pxsch.PRBSet) ~= length(usedsymbols)
            chname = extractBefore(ftable.MappingType,'Mapping');
            error('When PRBSet is a cell array, its length (%d) must be the same as the number of OFDM symbols defined for the %s allocation (%d).',...
                       length(pxsch.PRBSet),chname,length(usedsymbols));
        end
        allocatedPRB =  pxsch.PRBSet(usedsymbols);
    else
        allocatedPRB = pxsch.PRBSet;
    end
    
    % Capture the set of antenna ports associated with the DM-RS
    % Config type 1: p=(10)00-(10)07, Config type 2: p=(10)00-(10)11 (max ranges for PDSCH/PUSCH)
    if isfield(pxsch,'PortSet')
        ports = reshape(pxsch.PortSet,1,[]);
    else
        if isfield(pxsch,'NLayers') 
            ports = 0:pxsch.NLayers-1;     % Default to number of tx layers, if provided
        else
            ports = 0;                     % Otherwise single, single port
        end
        pxsch.PortSet = ports;
    end
    % Assign the PTRSPortSet with the lowest configured port, in the
    % absence of the field PTRSPortSet
    if ~isfield(pxsch,'PTRSPortSet') || isempty(pxsch.PTRSPortSet)
        pxsch.PTRSPortSet = min(ports);
    end

    % Assign the TransformPrecoding field to zero, in the absence of that
    % field
    if ~isfield(pxsch,'TransformPrecoding')
        pxsch.TransformPrecoding = 0;
    end

    % Capture the set of OFDM symbols (within the slot) that will carry the DM-RS
    % If not explicitly defined then calculate DM-RS containing OFDM symbols from parameter set
    if isfield(pxsch,'DMRSSymbolSet')
        dmrssymbols = intersect(pxsch.DMRSSymbolSet,allocatedsymbols);
        ldash = zeros(1,length(dmrssymbols));   % In the case, treat as all single symbols
    else
        % Establish whether it's type A or type B mapping
        mp = strcmpi(pxsch.(ftable.MappingType),["A","B"]);
        if ~any(mp)
            error('%s (%s) must be A or B.',ftable.MappingType,pxsch.(ftable.MappingType));
        end
        typeB = mp(2);
        [dmrssymbols,ldash] = getDMRSSymbolIndices(ftable.dmrsSymbolsTable,pxsch,typeB,allocatedsymbols);
    end
    
    % Detect and warn when no DM-RS symbols are defined for the input setup
    if ~isempty(allocatedsymbols) && isempty(dmrssymbols)
        chname = extractBefore(ftable.MappingType,'Mapping');
        warning('nr5g:DMRSParametersNoSymbols','No DM-RS symbols are defined for the given set of DM-RS configuration parameters and %s resource allocation. Cross-check with the %s DM-RS position tables.',...
                        chname,chname);
    end
    
    % Capture the slot number of the PDSCH/PUSCH
    % Used to locate the current slot in relation to any extended reserved resource patterns
    % and for the DM-RS sequence generation
    if ~isfield(pxsch,'NSlot')
        pxsch.NSlot = 0;     % Default to 0 if not present
    end

    % Capture the number of CDM groups without data
    if isfield(pxsch,'NumCDMGroupsWithoutData')
       cdmgroupsnodata = pxsch.NumCDMGroupsWithoutData;
       if ~any(cdmgroupsnodata == 0:3)
           error('NumCDMGroupsWithoutData (%d) must be one of (0,1,2,3).',cdmgroupsnodata);
       end
    else
       cdmgroupsnodata = 0;
    end
    
    % PRB level mapping (PRB carrying the PDSCH/PUSCH)
    % Find all allocated PRB (accounting for any reserved PRB) across
    % the allocated OFDM symbols
    
    % Create a cell array containing the PRB set used in each allocated OFDM symbol
    prbcell = cell(1,symbperslot);
    if iscell(allocatedPRB)
        prbcell(allocatedsymbols+1) = cellfun(@(x)reshape(mod(x,gnb.NRB),1,[]),allocatedPRB,'UniformOutput',false);
        prbset = prbcell{allocatedsymbols(1)+1}; % The initial value of shared channel allocation in the prbcell corresponds to PRBSet of first hop
    else
        prbset = reshape(mod(allocatedPRB,gnb.NRB),1,[]);        % Ensure that it is now a row for implicit expansion later
        prbcell(allocatedsymbols+1) = {prbset};
    end
       
    % Loop over the reserved resources and exclude them from the allocation PRB
    % across the allocated symbols, stored in prbcell
    if isfield(pxsch,'Reserved')
        for ri=1:length(pxsch.Reserved)

            % Reference the current reserved symbol/PRB indices pair and period
            reservedsymbols = pxsch.Reserved(ri).Symbols;
            reservedprb = reshape(pxsch.Reserved(ri).PRB,1,[]);
            reservedperiod = pxsch.Reserved(ri).Period;

            % Find any of the allocated symbols which overlap with reservations
            %
            % If the reserved period was empty then get number of complete slots
            % in reservation period and cyclically extend pattern to cover current slot 
            if isempty(reservedperiod)                
                reservedperiod = 0;
            end
            offset = mod(pxsch.NSlot,reservedperiod)*symbperslot;        % Symbol offset (whole number of slots) into pattern period                       
            inter = intersect(allocatedsymbols,reservedsymbols-offset);  % Find allocated symbols in current slot which contain reserved PRB
            % Reference the PRB associated with the overlapping symbols
            if iscell(reservedprb)
                if length(reservedprb) ~= length(reservedsymbols)  
                    error('Symbol-wise reserved PRB resources (defined in a cell array) and associated symbols should have the same number of elements.');
                end
                prbdiff = arrayfun(@(x)setdiff(prbcell{x+1},[reservedprb{x==(reservedsymbols-offset)}]), inter,'UniformOutput',false);
            else
                prbdiff = arrayfun(@(x)setdiff(prbcell{x+1},reservedprb), inter,'UniformOutput',false);
            end
            prbcell(inter+1) = prbdiff;
        end   
    end
    
    % VRB to PRB resource block mapping, TS 38.211 Section 7.3.1.6
    % Map the indices in prbcell if the PDSCH transmission requires
    % interleaved mapping 
    if isfield(pxsch,'VRBToPRBInterleaving') && pxsch.VRBToPRBInterleaving
        if ~isfield(pxsch,'VRBBundleSize') || isempty(pxsch.VRBBundleSize)
            L = 2; % Default
        else
            L = pxsch.VRBBundleSize;
        end
        if ~isfield(gnb,'RBOffset') || isempty(gnb.RBOffset)
            rboffset = 0; % Default
        else
            rboffset = gnb.RBOffset;
        end
        mapIndices = hPDSCHGenerateInterleavedIndices(gnb.NRB,rboffset,L);
        prbcell(~cellfun(@(x)isempty(x),prbcell)) = cellfun(@(x)mapIndices(ismember(mapIndices,x)),prbcell(~cellfun(@(x)isempty(x),prbcell)),'UniformOutput',false);
    end

    % RE level mapping (RE per PRB carrying the PDSCH/PUSCH)
    % Find all PDSCH/PUSCH RE per PRB for each OFDM symbol that will also carry DM-RS

    if ~any(pxsch.DMRSConfigurationType==[1,2])
        error('DMRSConfigurationType (%d) must be either 1 or 2.',pxsch.DMRSConfigurationType);
    end
    if pxsch.DMRSConfigurationType==1
        % Type 1: 6 DM-RS SC per PRB per CDM (every other SC)
        dmrssc = [0 2 4 6 8 10]';                   % RE indices in a PRB
        cdmgroups = mod(fix(ports/2),2);            % CDM groups associated with the used ports
        dshifts = cdmgroups;                        % Delta shift per used CDM group (0,1)... the *set* of shifts/groups that we will exclude
        dshiftsnodata = 0:min(cdmgroupsnodata,2)-1; % Delta shifts for CDM groups without data
    else      
        % Type 2: 4 DM-RS SC per PRB per CDM (2 groups of 2 SC)
        dmrssc = [0 1 6 7]';                            % RE indices in a PRB
        cdmgroups = mod(fix(ports/2),3);                % CDM groups associated with the used ports
        dshifts = 2*cdmgroups;                          % Delta shift per used CDM group (0,2,4)... the *set* of shifts/groups that we will exclude  
        dshiftsnodata = 2*(0:min(cdmgroupsnodata,3)-1); % Delta shifts for CDM groups without data
    end

    fullprb = ones(12,1);        % Binary map of all the subcarriers in a PRB
    dmrsre = dmrssc + [dshifts dshiftsnodata];% Implicit expansion of zero-shifted DM-RS RE positions (col) across all the shifts/CDMs (row) in play   
    fullprb(dmrsre+1) = 0;       % Clear all RE which will carry DM-RS in at least one port
    pxschre = find(fullprb)-1;   % Find PDSCH/PUSCH (non DM-RS) RE in a DM-RS containing symbol
    
    % Create the RE per PRB list across the allocated symbols, accounting
    % for DM-RS and non DM-RS carrying cases
    recell = cell(1,symbperslot);
    recell(allocatedsymbols+1) = {(0:11)'};  % Data RE per PRB in normal data symbols
    recell(dmrssymbols+1) = {pxschre};       % Data RE per PRB in DM-RS carrying symbols
      
    % Combine PRB oriented and RE per PRB oriented index arrays and expand (using implicit expansion)
    % into a column of linear indices for a single antenna/layer
    % Implicit expansion:                         Column             Row                                      Column
    slotindices = cell2mat( arrayfun(@(x) reshape(recell{x+1} + 12*prbcell{x+1} + 12*gnb.NRB*x,[],1), allocatedsymbols(:),'UniformOutput',false) );

    % Have the PT-RS indices and symbols through a local function provided
    % the PT-RS in enabled
    nslot = mod(pxsch.NSlot,slotsperframe);
    if ~isfield(pxsch,'EnablePTRS') || ~pxsch.EnablePTRS
        % PT-RS is disabled, provide the empty values
        ptrsSymbols = [];
        ptrsIndices = [];
        ptrsInfo = struct;
        ptrsInfo.PTRSSymbolSet = [];
    else
        % PT-RS is enabled
        if ~isfield(ftable,'getPTRSResources') || isempty(ftable.getPTRSResources)
            ftable.getPTRSResources = @ptrsCPOFDM;
        end
        [ptrsSymbols,ptrsIndices,ptrsInfo] = ptrsResources(ftable.getPTRSResources,pxsch,prbset,allocatedsymbols,dmrssymbols,nslot,gnb.NRB,symbperslot);
    end

    % Expand across all antenna ports/layers and make 1-based
    % Implicit expansion: Column                                        Row
    pxschIndices = slotindices(:) + (1 + (12*symbperslot*gnb.NRB)*(0:numel(ports)-1));

    % For CP-OFDM, remove the ptrsIndices from the slotIndices thereby
    % reducing the symbol capacity. For DFT-s-OFDM, use the slotIndices
    % without removing the ptrsIndices, since, both data and PT-RS are
    % mapped within the shared channel indices itself
    if ~pxsch.TransformPrecoding
        [~,ptrsIndGrid] = nrExtractResources(ptrsIndices,zeros([gnb.NRB*12 symbperslot numel(ports)]));
        pxschIndices = reshape(setdiff(pxschIndices,ptrsIndGrid,'stable'),[],numel(ports));

        % Remove PT-RS indices and symbols from the reserved resource
        % blocks
        [~,validPTRS] = intersect(ptrsIndGrid(:,1),slotindices+1,'stable');
        ptrsIndices = reshape(ptrsIndices(validPTRS),[],size(ptrsIndices,2));
        ptrsSymbols = reshape(ptrsSymbols(validPTRS),[],size(ptrsSymbols,2));
    end

    % Channel bit capacity information
    nports = numel(ports);                      % Numbers of ports/layers across all codewords
    ncw = 1 + (nports > 4);                     % Number of codewords, deduced from total layers
    nlayers = fix((nports + (0:ncw-1))/ncw);    % Number of layers per codeword
    pxschInfo.Gd = size(pxschIndices,1);        % Number of QAM data symbols in one layer/antenna grid
    % For DFT-s-OFDM, the PUSCH symbol capacity must be reduced by PT-RS
    % symbol capacity, rather than removing the PT-RS indices from
    % slotIndices
    if pxsch.TransformPrecoding
        pxschInfo.Gd = size(pxschIndices,1) - size(ptrsIndices,1); % Number of QAM data symbols in one layer/antenna grid
    end

    % Validate modulation type and translate into bits per symbol
    fullmodlist = ["pi/2-BPSK","BPSK","QPSK","16QAM","64QAM","256QAM"]'; % All NR mod schemes
    modulation = reshape(string(pxsch.Modulation),1,[]);              % Turn into a string row
    qm = [1 1 2 4 6 8]*(lower(modulation) == lower(fullmodlist));     % Test each element against all mods then multiply against bps
    if ~all(qm)
        error("The modulation (%s) must be one of the set (%s).",join(modulation(logical(qm == 0))), join(fullmodlist',','));
    end
    pxschInfo.G = qm.*pxschInfo.Gd.*nlayers;     % And scale by layers/ports and modulation scheme
    pxschInfo.NREPerPRB  = 12*(length(allocatedsymbols)-length(dmrssymbols)) + length(pxschre)*length(dmrssymbols);
  
    % Additional DM-RS related information
    pxschInfo.DMRSSymbolSet = dmrssymbols;      % OFDM symbols containing DM-RS
    pxschInfo.CDMGroups = cdmgroups;            % CDM group associated with each port

    % Additional PT-RS related information
    pxschInfo.PTRSSymbolSet = ptrsInfo.PTRSSymbolSet;

    % Create the DM-RS QPSK symbols and resource element indices associated
    % with the shared transmission
    %
    % First create *unshifted* version of the DM-RS RE indices in a single antenna plane
    % Implicit expansion:                           Column        Row          Scalar offset          Column
    dmslotindices = cell2mat( arrayfun(@(x) reshape(dmrssc + 12*prbcell{x+1} + 12*gnb.NRB*x,[],1), dmrssymbols(:),...
                              'UniformOutput',false) );
    
    % Expand across all antenna ports/layers, applying DM-RS shifts for each port, and make 1-based
    % Implicit expansion: Column             Row          Scalar Offset          Row
    dmrsIndices = dmslotindices(:) + (1 + dshifts + (12*symbperslot*gnb.NRB)*(0:nports-1));
    
    % If an overriding DM-RS sequence generator is defined then use it
    % otherwise use the local NIDNSCID/NSCID PRBS based generator
    if ~isfield(ftable,'getDMRSSequence') || isempty(ftable.getDMRSSequence)
        ftable.getDMRSSequence = @prbsDMRSSequence;
    end
    % Create the matching set of DM-RS QPSK symbols for all ports, and 
    % record the frequency and time cover codes used for each port
    % First calculate the frame relative slot number from the (absolute) input slot number
    [dmrSymbols,covers] = getDMRSValues(ftable.getDMRSSequence,pxsch,prbcell(dmrssymbols+1),nslot,dmrssymbols,ldash,symbperslot);
    
    % Calculate CDM despreading lengths (for use in nrChannelEstimate)
    pxschInfo.CDMLengths = getCDMLengths(cdmgroups,covers);
    
    % Additional structural information to be passed back to the calling function
    genInfo.SymbolsPerSlot = symbperslot;
    genInfo.PortSet = ports;
end

% Get the OFDM symbol indices containing DM-RS for the PDSCH/PUSCH allocation
% Uses a function handle callback to look up the initial l^bar position values
% which are then adjusted for the mapping type and duration (single/double)
function [dmrssymbols,ldash] = getDMRSSymbolIndices(gettablesymbols,pxsch,typeb,allocatedsymbols)  
    
    % Check DMRSLength value for later use
    dmslength = pxsch.DMRSLength;
    if ~isscalar(dmslength) || ~any(dmslength == [1 2])
        error('DMRSLength (%s) must be either 1 or 2.', join(string(dmslength),','));
    end

    % Get PDSCH/PUSCH allocation duration, l_d  
    [lb,ub] = bounds(allocatedsymbols);
    if ~typeb
        lb = 0;
    end
    nsymbols = ub - lb + 1;
    
    if nsymbols
        % Get the OFDM symbol indices, l^bar, that will carry DM-RS for the channel instance
        dmrssymbols = gettablesymbols(pxsch,typeb,nsymbols);
    else
        dmrssymbols = [];
    end
    
    % l' values associated with DM-RS symbol indices
    ldash = zeros(1,length(dmrssymbols)*dmslength);
    % Adjust table information
    if ~isempty(dmrssymbols)
        % Adjust indices for the relative offset of the mapping type
        if typeb
           dmrssymbols = dmrssymbols+allocatedsymbols(1);   % If type B (non slot-wise)
        else
           dmrssymbols(1) = pxsch.DMRSTypeAPosition;        % If type A (slot-wise) then 2 or 3
        end
        % Adjust for double-symbol DM-RS
        % In the operational system, if RRC configured with max length 2
        % then the actual length is DCI signaled
        if pxsch.DMRSLength == 2
            dmrssymbols = reshape([dmrssymbols; dmrssymbols+1],1,[]);
            ldash(2:2:end) = 1;
        end
        % For non-standard set-ups, only return the DM-RS symbol indices that 
        % overlap the actual allocation indices
        [dmrssymbols,~,indices] = intersect(allocatedsymbols,dmrssymbols);
        ldash = ldash(indices);
    end
end

% Get all complex symbol values associated with the DM-RS for the transmission
% Uses a function handle callback to 
% Each separate antenna port is a column of the returned matrix which concatenates
% all the DM-RS values for all DM-RS carrying OFDM symbols associated with transmission
function [symb,covers] = getDMRSValues(getDMRSSequence,pxsch,prbset,nslot,dmrssymbols,ldashvalues,symbperslot)
       
    % Capture the PRB/subcarrier reference point associated with PRB set values, 
    % for the DM-RS sequence indexing
    if isfield(pxsch,'PRBRefPoint')
       prbrefpoint = pxsch.PRBRefPoint;
    else
       prbrefpoint = 0;
    end
    
    % Capture set of antenna ports required
    ports = pxsch.PortSet;
    
    % 6 DM-RS QPSK symbols (type 1) or 4 DM-RS QPSK symbols (type 2) per PRB
    ndmrsre = 4 + (pxsch.DMRSConfigurationType==1)*2;     
    
    % Loop over the OFDM symbols containing DM-RS 
    symcell = cell(1,length(dmrssymbols));
    for i=1:length(dmrssymbols)
        symcell{i} = reshape(getDMRSSequence(pxsch,ndmrsre,prbset{i},prbrefpoint,nslot,dmrssymbols(i),ldashvalues(i),symbperslot),1,[]);
    end

    % Accumulated lengths of the DM-RS when concatenated across the OFDM symbols 
    cndmrs = [0 cumsum(cellfun(@length,symcell))];
    % Preallocate array for the returned DM-RS
    symb = zeros(cndmrs(end),numel(ports));

    % Establish the max number of DM-RS ports, depending on the config type
    maxnports = 8+(pxsch.DMRSConfigurationType==2)*4;
    % Loop over the ports and apply masks, recording covers used in mask
    % values
    covers = ones(numel(ports),2);
    for pidx = 1:numel(ports)        
        % Mask applied if (p is odd (f part)) or (p >= (100)4 (type 1) or (100)6 (type 2) and double-symbol (t part))
        % Test for 'f' part of mask (subcarrier, frequency part)
        pv = ports(pidx);
        if mod(pv,2)
          fmask = [1 -1];
          covers(pidx,1) = -1;
        else
          fmask = 1;
        end
        % Test for the 't' part, based on the l' values
        tmask = 1-2*(ldashvalues*(pv >= maxnports/2));
        if (any(tmask==-1))
            covers(pidx,2) = -1;
        end
        % Apply combined time and frequency mask values and concatenate 
        % across all the OFDM symbols
        for i = 1:length(symcell)             
            symb(cndmrs(i)+1:cndmrs(i+1),pidx) = symcell{i}.*repmat(tmask(i)*fmask,size(symcell{i})./size(fmask));
        end
    end
end

% Generate PRBS based DM-RS sequence for a single OFDM symbol
function symbols = prbsDMRSSequence(pxsch,ndmrssc,prbset,prbrefpoint,nslot,nsymbol,ldash,symbperslot) %#ok<INUSL>
    
    % Cache the scrambling IDs
    nidnscid = pxsch.NIDNSCID;
    nscid = pxsch.NSCID;
    
    if ~isempty(prbset)
        % Generate PRBS for DM-RS sequence which covers the PRB allocation set range 
        [minprb,maxprb] = bounds(prbset);
        cinit = mod(2^17*(symbperslot*nslot + nsymbol + 1)*(2*nidnscid + 1) + 2*nidnscid + nscid,2^31);
        prbs = reshape(nrPRBS(cinit,2*ndmrssc*[prbrefpoint+minprb maxprb-minprb+1]),2*ndmrssc,[]);
    
        % Extract PRBS values associated with PRB and turn into complex DM-RS symbols
        bpsk = 1/sqrt(2)*(1-2*reshape(prbs(:,prbset-minprb+1),2,[])');
        symbols = complex(bpsk(:,1),bpsk(:,2));
    else
        symbols = complex([]);
    end
end   

% Generate 0-based PRB indices for a PDSCH transmission. PDSCH
% transmissions not scheduled with DCI format 1-0, as specified in
% TS 38.211 Section 7.3.1.6.
function PRBindices = hPDSCHGenerateInterleavedIndices(nrb,rboffset,L)
    
	rboffsetModL = mod(rboffset,L);
	
    % PRB bundle generation
    Nb = ceil((nrb+rboffsetModL)/L); % Number of PRB bundles
    
    % If only one bundle, force the indices to be [0:NRB-1]
    if Nb == 1
        PRBindices = (0:nrb-1);
        return
    end
    
    % Generate the 0-based RB bundle indices
    numRBinBundle = zeros(1,Nb);
    numRBinBundle(1) = L-rboffsetModL;
    numRBinBundle(end) = mod(rboffset+nrb,L);
    if ~numRBinBundle(end)
        numRBinBundle(end) = L;
    end
    numRBinBundle(2:end-1) = L;
    
    % Generate the 0-based PRB bundles indices for the interleaved mapping
    R = 2;
    C = floor(Nb/R);
    r = (0:R-1)';
    c = 0:C-1;
    prbbInd = r*C+c;    % Indices of the PRB bundle
    prbbInd = prbbInd(:)';
    prbbInd(Nb) = Nb-1; % Last VRB bundle is mapped to the last PRB bundle
    
    % Note that prbbInd(1)=0 and prbbInd(end)=Nb-1 always by design.
    % Thus, only prbbInd(2:end-1) are actively considered.
	prbInd = (prbbInd(2:end-1).*numRBinBundle(2:end-1)) - rboffsetModL;
    
    % Extend the bundle indices to the actual RB indices
    PRBindices = NaN(1,nrb);
    PRBindices(1:numRBinBundle(1)) = (0:numRBinBundle(1)-1); % First bundle
    PRBindices_tmp = prbInd+(0:L-1)';
    PRBindices(1+numRBinBundle(1):end-numRBinBundle(end)) = PRBindices_tmp(:)'; % Bundles in the middle
    PRBindices(end-numRBinBundle(end)+1:end) = max(PRBindices)+1+(0:numRBinBundle(end)-1); % Last bundle
    
end

% Calculate CDM despreading lengths (for use in nrChannelEstimate), by
% finding the maximum number of different cover codes used across time
% and frequency for each CDM group
function cdmLengths = getCDMLengths(cdmgroups,covers)
    
    ug = unique(cdmgroups);
    
    cdmLengths = zeros(numel(ug),2);
    
    for gi = 1:numel(ug)
        ci = find(cdmgroups==ug(gi));
        cdmLengths(gi,1) = numel(unique(covers(ci,1)));
        cdmLengths(gi,2) = numel(unique(covers(ci,2)));
    end
    
    cdmLengths = max(cdmLengths,[],1);
     
end

function [ptrsSym,ptrsInd,info] = ptrsResources(getPTRSResources,pxsch,prbSet,allocatedSymbols,dmrsSymbols,nslot,nsizegrid,symbperslot)
%ptrsResources Provides the PT-RS resources based on the waveform
% along with structural information

    % Initialize the pxsch fields that affect the shared channel processing
    if ~isfield(pxsch,'RNTI')
        pxsch.RNTI = 0;
    end
    if ~isfield(pxsch,'IntraSlotFreqHopping')
        pxsch.IntraSlotFreqHopping = 'disabled';
    end

    % Symbol allocation of shared channel
    symbolSet = unique(allocatedSymbols);
    pxschSymAllocation = [symbolSet(1) length(symbolSet)];

    % Assign time density (lptrs) with the default value, if not provided
    if ~isfield(pxsch,'PTRSTimeDensity')
        pxsch.PTRSTimeDensity = 1;
    end
    lptrs = double(pxsch.PTRSTimeDensity);

    % Get PT-RS symbol indices
    if ~isempty(dmrsSymbols) && ~isempty(lptrs)
        ptrsSymInd = pxschSymAllocation(1):lptrs:dmrsSymbols(1)-1;
        for i = 1:numel(dmrsSymbols)
            ptrsSymInd(ptrsSymInd >= dmrsSymbols(i)) = [];
            ptrsSymInd = [ptrsSymInd (dmrsSymbols(i)+lptrs):lptrs:(sum(pxschSymAllocation)-1)]; %#ok<AGROW>
        end
    else
        ptrsSymInd = zeros(1,0);
    end

    % Initialize outputs and intermediate variables
    ptrsSym = complex(zeros(0,1));
    ptrsInd = zeros(0,1);

    % Assign frequency density, number of PT-RS samples, and number of
    % PT-RS groups with default values, if the field(s) are not provided
    nRB = numel(unique(prbSet));
    if pxsch.TransformPrecoding
        if ~isfield(pxsch,'PTRSNumSamples')
            pxsch.PTRSNumSamples = 2;
        end
        if ~isfield(pxsch,'PTRSNumGroups')
            pxsch.PTRSNumGroups = 2;
        end
        % Don't generate PT-RS when number of subcarriers is less than the
        % number of PT-RS symbols and when any of PTRSNumSamples and
        % PTRSNumGroups is empty
        cond = isempty(pxsch.PTRSNumSamples) || isempty(pxsch.PTRSNumGroups) ...
            || ((nRB*12)<(pxsch.PTRSNumSamples*pxsch.PTRSNumGroups));
        % Assign the PTRSNumSamples and PTRSNumGroups to freqDensity to use
        % it directly in getPTRSResources
        freqDensity = [pxsch.PTRSNumSamples pxsch.PTRSNumGroups];
    else
        if ~isfield(pxsch,'PTRSFrequencyDensity')
            pxsch.PTRSFrequencyDensity = 2;
        end
        if ~isfield(pxsch,'PTRSREOffset')
            pxsch.PTRSREOffset = '00';
        end
        % Don't generate PT-RS when FrequencyDensity is empty
        cond = isempty(pxsch.PTRSFrequencyDensity);
        % Assign the PTRSFrequencyDensity to freqDensity to use it directly in
        % getPTRSResources
        freqDensity = pxsch.PTRSFrequencyDensity;
    end

    if ~(isempty(prbSet) || cond || ~numel(ptrsSymInd))
        % Validate PT-RS configuration
        validatePTRS(pxsch);

        % Get the PT-RS symbols and indices
        [ptrsSym,ptrsInd] = getPTRSResources(pxsch,prbSet,pxschSymAllocation,dmrsSymbols,ptrsSymInd,freqDensity,nslot,symbperslot,nsizegrid);
    end

    % Output structure info
    info = struct;
    info.DMRSSymbolSet = dmrsSymbols;
    info.PTRSSymbolSet = ptrsSymInd;

end

function validatePTRS(pxsch)
% Validate the PT-RS configuration parameters and dependencies

    % Check the time density
    if ~any(pxsch.PTRSTimeDensity == [1 2 4])
        error('nr5g:PTRS:InvalidTimeDensity','The time density of PT-RS (%d) must be 1, 2, or 4.',pxsch.PTRSTimeDensity);
    end

    if ~pxsch.TransformPrecoding
        % Check the frequency density
        if ~any(pxsch.PTRSFrequencyDensity == [2 4])
            error('nr5g:PTRS:InvalidFrequencyDensity','The frequency density of PT-RS (%d) must be either 2 or 4.',pxsch.PTRSFrequencyDensity);
        end

        % Check the PTRSREOffset
        validatestring(pxsch.PTRSREOffset,{'00','01','10','11'},'','REOffset of PT-RS');

        % Check if DM-RS antenna ports are in valid range while
        % transmitting PT-RS
        maxDMRSPort = max(pxsch.PortSet);
        maxAcceptable = 2*pxsch.DMRSConfigurationType+1;
        if maxDMRSPort > maxAcceptable
            error('nr5g:PTRS:InvalidDMRSPorts','For DM-RS configuration type %d, the maximum DM-RS antenna port (%d) exceeds the possible range 0 to %d, when PT-RS is configured.',pxsch.DMRSConfigurationType,maxDMRSPort,maxAcceptable)
        end

        ptrsPorts = unique(reshape(pxsch.PTRSPortSet,1,[]));
        ptrsPortsLen = length(ptrsPorts);
        if ~(ptrsPortsLen == 1 || ptrsPortsLen == 2)
            error('nr5g:PTRS:InvalidNumPTRSPorts','The number of unique PT-RS ports (%d) must be 1 or 2.',ptrsPortsLen);
        end

        % Check if PT-RS port set is a subset of DM-RS port set
        if ptrsPortsLen == 2
            if ~any(ptrsPorts(1) == pxsch.PortSet) || ~any(ptrsPorts(2) == pxsch.PortSet)
                error('nr5g:PTRS:InvalidPTRSPorts','The PTRSPortSet [%d,%d] must be a subset of PortSet.',...
                    ptrsPorts(1), ptrsPorts(2));
            end
        else
            if ~any(ptrsPorts(1) == pxsch.PortSet)
                error('nr5g:PTRS:InvalidPTRSPort','The PTRSPortSet (%d) must be a subset of PortSet.',ptrsPorts(1));
            end
        end
    else
        % Check the PTRSNumSamples and PTRSNumGroups for DFT-s-OFDM
        if ~any(pxsch.PTRSNumSamples == [2 4])
            error('nr5g:PTRS:InvalidNumSamples','The number of PT-RS samples (%d) must be either 2 or 4.',pxsch.PTRSNumSamples);
        end
        if ~any(pxsch.PTRSNumGroups == [2 4 8])
            error('nr5g:PTRS:InvalidNumGroups','The number of PT-RS groups (%d) must be 2, 4, or 8.',pxsch.PTRSNumGroups);
        end
        if pxsch.PTRSNumSamples == 2 && pxsch.PTRSNumGroups == 8
            error('nr5g:PTRS:InvalidCombination','When the number of PT-RS groups is 8, the number of PT-RS samples (2) must be 4.');
        end
    end

end

function [kRefTable,dmrsSCPattern,nDMRSSC] = ptrsSubCarrierInfo(dmrsConfigType)
%ptrsSubcarrierInfo Provides the subcarrier related information for PT-RS
% processing

    if dmrsConfigType == 1
        % Table 7.4.1.2.2-1/6.4.1.2.2.1-1, TS 38.211
        kRefTable = [0 2 6 8;...
                     2 4 8 10;...
                     1 3 7 9;...
                     3 5 9 11];
        % DM-RS subcarrier locations pattern per RB per port
        dmrsSCPattern = [0:2:10; 1:2:11];
        % Number of subcarriers carrying DM-RS symbols per RB per port
        nDMRSSC   = 6;
    else
        % Table 7.4.1.2.2-1/6.4.1.2.2.1-1, TS 38.211
        kRefTable = [0 1 6 7;...
                     1 6 7 0;...
                     2 3 8 9;...
                     3 8 9 2;...
                     4 5 10 11;...
                     5 10 11 4];
        % DM-RS subcarrier locations pattern per RB per port
        dmrsSCPattern = [0 1 6 7; 2 3 8 9; 4 5 10 11];
        % Number of subcarriers carrying DM-RS symbols per RB per port
        nDMRSSC   = 4;
    end
end

function [sym,ind] = ptrsCPOFDM(pxsch,prbSet,pxschSymAllocation,dmrsSymbols,ptrsSymInd,kPTRS,nslot,symbperslot,nsizegrid)
%ptrsCPOFDM Provides PT-RS symbols and indices for CP-OFDM waveform

    % PRB reference point to account for the DM-RS symbols from that point
    if isfield(pxsch,'PRBRefPoint')
       prbrefpoint = pxsch.PRBRefPoint;
    else
       prbrefpoint = 0;
    end

    % Assign the DM-RS port set to a variable
    dmrsPortSet = double(pxsch.PortSet);

    % Number of subcarriers
    nRBSC      = 12;              % Number of subcarrier in a resource block
    nCarrierSC = nRBSC*nsizegrid; % Number of subcarriers in a carrier

    % Number of OFDM symbols allocated for PXSCH
    nPXSCHSymbols = pxschSymAllocation(2);

    % Subcarriers related information
    [kRefTable,dmrsSCPattern,nDMRSSC] = ptrsSubCarrierInfo(pxsch.DMRSConfigurationType);

    % Resource element offset kRERef of each port
    ptrsPorts = unique(pxsch.PTRSPortSet);
    colIndex = strcmpi(pxsch.PTRSREOffset,{'00','01','10','11'});
    kRERef = kRefTable(ptrsPorts+1,colIndex);

    % Initialize the variables
    nPTRSPorts = numel(ptrsPorts);
    ptrsTemp = cell(1,nPTRSPorts);
    indTemp = cell(1,nPTRSPorts); % Column cell type

    % Number of PT-RS ports
    numPTRSPortsEquals2 = (nPTRSPorts==2);

    % PRB set of PXSCH allocation
    prbset = unique(reshape(prbSet,1,[]));
    nPXSCHRB = numel(prbset);
    nPXSCHSC = nPXSCHRB*12;

    % DM-RS symbol indices in second hop, intra-slot frequency hopping flag
    dmrsHop2 = dmrsSymbols(dmrsSymbols(1,:) >= (floor(nPXSCHSymbols/2)+pxschSymAllocation(1)));
    dmrsHop2Flag = any(dmrsHop2);
    intraSlotFreqEnabled = strcmpi(pxsch.IntraSlotFreqHopping,'enabled');

    % PRB set in each hop
    prbsetHop = zeros(1+(dmrsHop2Flag && intraSlotFreqEnabled),nPXSCHRB); % Each row corresponds to a hop
    prbsetHop(1,:) = prbset;

    % DM-RS symbols in each hop
    dmrsSym = complex(zeros(nDMRSSC,nPXSCHRB,1+(dmrsHop2Flag && intraSlotFreqEnabled)));
    dmrsSymTemp = prbsDMRSSequence(struct('NIDNSCID',pxsch.NIDNSCID,'NSCID',pxsch.NSCID),...
        nDMRSSC,prbset,prbrefpoint,nslot,dmrsSymbols(1),0,symbperslot);
    dmrsSym(:,:,1) = reshape(dmrsSymTemp,[],nPXSCHRB);

    if dmrsHop2Flag && intraSlotFreqEnabled
        % Provide DM-RS symbols only when intra-slot frequency hopping is
        % enabled and at-least one DM-RS symbol index is present in second hop
        prbsetHop(2,:) = prbset-min(prbset)+pxsch.RBOffset;
        dmrsSymTemp = prbsDMRSSequence(struct('NIDNSCID',pxsch.NIDNSCID,'NSCID',pxsch.NSCID),...
            nDMRSSC,prbsetHop(2,:),prbrefpoint,nslot,dmrsHop2(1),0,symbperslot);
        dmrsSym(:,:,2) = reshape(dmrsSymTemp,[],nPXSCHRB);
    end

    % Resource block offset kRBRef
    if mod(nPXSCHRB,kPTRS) == 0
        kRBRef = mod(pxsch.RNTI,kPTRS);
    else
        kRBRef = mod(pxsch.RNTI,mod(nPXSCHRB,kPTRS));
    end

    for p = 1:nPTRSPorts
        % PT-RS symbols for each hop
        [~,scIndex] = find(repmat(kRERef(p),size(dmrsSCPattern)) == dmrsSCPattern);
        ptrsHop = permute(dmrsSym(scIndex,kRBRef(1)+1:kPTRS:end,:),[2 3 1]); % Each column corresponds to a hop

        % PT-RS indices in second hop
        SecondHopInd = ptrsSymInd >= floor(nPXSCHSymbols/2);

        % PT-RS symbols and indices for each port
        ip = 0:floor((nPXSCHRB-((1+kRERef(p))/nPXSCHSC)-kRBRef(1))/kPTRS);
        if intraSlotFreqEnabled
            % Intra-slot frequency hopping enabled
            ptrsTemp{1,p} = repmat(ptrsHop(:,1),nnz(~SecondHopInd),1);
            symbolOffsetFirstHop = nCarrierSC*(ptrsSymInd(~SecondHopInd));
            if ~isempty(symbolOffsetFirstHop)
                indFirstHop = reshape(reshape(kRERef(p)+(prbsetHop(1,ip*kPTRS+1)+kRBRef(1))*nRBSC,[],1) + symbolOffsetFirstHop, [], 1);
                indTemp{1,p} = indFirstHop;
            end
            if dmrsHop2Flag
                % Provide PT-RS symbols and indices for second hop, when
                % at-least one DM-RS symbol index is present in second hop
                ptrsTemp{1,p} = [ptrsTemp{1,p}; repmat(ptrsHop(:,2),nnz(SecondHopInd),1)];
                symbolOffsetSecondHop = nCarrierSC*(ptrsSymInd(SecondHopInd));
                if ~isempty(symbolOffsetSecondHop)
                    indSecondHop = reshape(reshape(kRERef(p)+(prbsetHop(2,ip*kPTRS+1)+kRBRef(1))*nRBSC,[],1) + symbolOffsetSecondHop, [], 1);
                    indTemp{1,p} = [indTemp{1,p}; indSecondHop];
                end
            end
        else
            % Intra-slot frequency hopping disabled
            ptrsTemp{1,p} = repmat(ptrsHop,numel(ptrsSymInd),1);
            symbolOffsets = nCarrierSC*ptrsSymInd;
            indTemp{1,p} = reshape(reshape(kRERef(p)+(prbsetHop(1,ip*kPTRS+1)+kRBRef(1))*nRBSC,[],1) + symbolOffsets, [], 1);
        end
    end

    % NonOrthogonal port multiplexing
    [~,Offset1] = max(dmrsPortSet == ptrsPorts(1)); % Layer corresponding to association of PT-RS port 1 with DM-RS port
    if numPTRSPortsEquals2
        % Layer corresponding to association of PT-RS port 2 with DM-RS port
        [~,Offset2] = max(dmrsPortSet == ptrsPorts(2));
    else
        Offset2 = [];
    end
    sym = cell2mat(ptrsTemp);
    Offset = [Offset1 Offset2];
    ind = 1+cell2mat(indTemp)+nCarrierSC*symbperslot*(Offset-1);

end
