classdef mynrCarrierConfig < nrCarrierConfig
    %mynrCarrierConfig nrCarrierConfig with a custom PRB width.
    %
    %   The class keeps all nrCarrierConfig behavior and adds only one
    %   property, SubcarriersPerPRB. Standard toolbox functions are still
    %   used whenever their inputs do not assume 12 subcarriers per PRB.

    properties
        SubcarriersPerPRB (1,1) {mustBeInteger,mustBePositive} = 12
    end

    properties (Dependent)
        NCustomSubcarriers
    end

    methods
        function obj = mynrCarrierConfig(varargin)
            obj@nrCarrierConfig;
            if mod(numel(varargin),2) ~= 0
                error('mynrCarrierConfig:InvalidNameValue','Name-value inputs must appear in pairs.');
            end
            for i = 1:2:numel(varargin)
                obj.(char(varargin{i})) = varargin{i+1};
            end
        end

        function k = get.NCustomSubcarriers(obj)
            k = obj.NSizeGrid * obj.SubcarriersPerPRB;
        end
    end
end
