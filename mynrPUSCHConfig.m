classdef mynrPUSCHConfig < nr5g.internal.BWPSizeStart & myConfigBase %nr5g.internal.pusch.ConfigBase
    %nrPUSCHConfig PUSCH configuration object
    %   The class keeps all nr5g.internal.BWPSizeStart behavior except the
    %   nr5g.internal.pusch.ConfigBase which includes the DMRS
    %   configuration. Standard toolbox functions are still used whenever
    %   their inputs do not assume 12 subcarriers per PRB.
    %
    %
    %   nrPUSCHConfig properties (configurable):
    %
    %   NSizeBWP           - Size of the bandwidth part (BWP) in terms
    %                        of number of physical resource blocks
    %                        (PRBs) (1...275) (default [])
    %   NStartBWP          - Starting PRB index of BWP relative to
    %                        common resource block 0 (CRB 0) (0...2473)
    %                        (default [])
    %   Modulation         - Modulation scheme(s) of codeword(s)
    %                        ('QPSK' (default), 'pi/2-BPSK', '16QAM', '64QAM', '256QAM')
    %   NumLayers          - Number of transmission layers (1...8) (default 1)
    %   MappingType        - PUSCH mapping type ('A' (default), 'B')
    %   SymbolAllocation   - OFDM symbol allocation of PUSCH within a slot
    %                        (default [0 14])
    %   PRBSet             - PRBs allocated for PUSCH within the BWP (default 0:51)
    %   TransformPrecoding - Flag to enable transform precoding (0(default), 1).
    %                        0 indicates that transform precoding is
    %                        disabled and the waveform type is CP-OFDM. 1
    %                        indicates that transform precoding is enabled
    %                        and the waveform type is DFT-s-OFDM
    %   TransmissionScheme - PUSCH transmission scheme ('nonCodebook' (default), 'codebook')
    %   NumAntennaPorts    - Number of antenna ports (1 (default), 2, 4, 8)
    %   TPMI               - Transmitted precoding matrix indicator (0...304) (default 0)
    %   CodebookType       - Codebook type ('codebook1_ng1n4n1' (default),
    %                        'codebook1_ng1n2n2', 'codebook2', 'codebook3', 'codebook4')
    %   FrequencyHopping   - Flag to enable frequency hopping
    %                        ('neither' (default), 'intraSlot', 'interSlot')
    %   SecondHopStartPRB  - Starting PRB of second hop relative to the
    %                        BWP (0...274) (default 1)
    %   Interlacing        - Enable interlacing (default false)
    %   RBSetIndex         - Resource block set index (default 0)
    %   InterlaceIndex     - Interlace indices (0...9) (default 0)
    %   BetaOffsetACK      - Beta offset for HARQ-ACK (default 20)
    %   BetaOffsetCSI1     - Beta offset for CSI part 1 (default 6.25)
    %   BetaOffsetCSI2     - Beta offset for CSI part 2 (default 6.25)
    %   UCIScaling         - Scaling factor to limit the number of resource
    %                        elements for UCI on PUSCH (default 1)
    %   NID                - PUSCH scrambling identity (0...1023) (default [])
    %   RNTI               - Radio network temporary identifier (0...65535)
    %                        (default 1)
    %   NRAPID             - Random access preamble index to initialize the
    %                        scrambling sequence for msgA on PUSCH (0...63)
    %                        (default [])
    %   DMRS               - PUSCH-specific DM-RS configuration object
    %   EnablePTRS         - Enable or disable the PT-RS configuration
    %                        (0 (default), 1)
    %   PTRS               - PUSCH-specific PT-RS configuration object
    %
    %   nrPUSCHConfig properties (read-only):
    %
    %   NumCodewords          - Number of codewords
    % 
    %   nrPUSCHConfig methods:
    %   
    %   nrTBS - Transport block size(s) associated with transmission

    % Public, tunable properties
    properties
        %NID Physical shared channel scrambling identity
        % Specify the physical shared channel scrambling identity as a
        % scalar nonnegative integer. The value must be in the range
        % 0...1023. It is the dataScramblingIdentityPUSCH (0...1023), if
        % configured, else it is the physical layer cell identity
        % (0...1007). Use empty ([]) to make this property equal to the
        % <a href="matlab:help('nrCarrierConfig/NCellID')"
        % >NCellID</a> property of nrCarrierConfig. The default value is [].
        NID = [];
    end
    
    properties (Hidden)
        % Custom property list to change the order of display properties
        CustomPropList = {'NSizeBWP','NStartBWP','Modulation','NumLayers',...
            'MappingType','SymbolAllocation','PRBSet',...
            'TransformPrecoding','TransmissionScheme','NumAntennaPorts',...
            'TPMI','CodebookType','FrequencyHopping','SecondHopStartPRB',...
            'Interlacing','RBSetIndex','InterlaceIndex',...
            'BetaOffsetACK','BetaOffsetCSI1','BetaOffsetCSI2','UCIScaling',...
            'NID','RNTI','NRAPID','DMRS','EnablePTRS','PTRS','NumCodewords'};
    end

    methods
        % Constructor
        function obj = mynrPUSCHConfig(varargin)
            % Get the value of NStartBWP from the name-value pairs
            nStartBWP = nr5g.internal.parseProp('NStartBWP',[],varargin{:});
            % Get the value of NSizeBWP from the name-value pairs
            nSizeBWP = nr5g.internal.parseProp('NSizeBWP',[],varargin{:});
            % Get the value of NID from the name-value pairs
            nid = nr5g.internal.parseProp('NID',[],varargin{:});
            % Get the value of NRAPID from the name-value pairs
            nrapid = nr5g.internal.parseProp('NRAPID',[],varargin{:});
            
            % Support name-value pair arguments when constructing object
            obj@myConfigBase(...
                'NStartBWP', nStartBWP, ...
                'NSizeBWP', nSizeBWP, ...
                'NID', nid, ...
                'NRAPID', nrapid, ...
                varargin{:});
        end
        
        % Self-validate and set properties
        function obj = set.NID(obj,val)
            prop = 'NID';
            nid = validateNID(obj,val);
            obj.(prop) = nid;
        end

        function tbs = nrTBS(obj,tcr,varargin)
            % Signatures supported are,
            % nrTBS(pusch, tcr)
            % nrTBS(pusch, tcr,xOh)
            % nrTBS(pusch, tcr,xOh, slotcarrier) - Mandatory in the interlaced PRB allocation case

            % Input argument checking and preprocessing 
            narginchk(1,4);
            if nargin == 4
                validateattributes(varargin{2},"nrCarrierConfig","scalar","nrTBS","CARRIER");
            end

            [tcr,xOh] = nr5g.internal.TBSDetermination.transportChParamsLinear(tcr, 0, varargin{:});  % Don't include TB scaling in the argument parsing  
            
            interlacing = obj.Interlacing;
            if interlacing
                % Check input signature for nrTBS(channel,tcr,xoh,carrier:nrCarrierConfig) syntax
                errFlag = interlacing && nargin ~= 4;
                coder.internal.errorIf(errFlag,'nr5g:nrTBS:InvalidSigForInterlacedPUSCH');
    
                % Dispatch to interlaced PUSCH specific function to include the PRB
                % allocation for this case
                tbs = getTBSEntryInterlaced(obj,tcr,xOh,obj,varargin{2});
            else            
                % Otherwise perform standard calculation
                tbs = nr5g.internal.TBSDetermination.getTBSEntry(obj,tcr,xOh,1);
            end
        end
       
    end

end
