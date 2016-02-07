function [int, varargout] = get_section_id(sect, varargin)
%GET_SECTION_ID find the corresponding section ID by name
%
%
%Marcell Stippinger, 2016

nXtra = nargin - 1;

if strcmp(sect,'mid_arm')
    int = 1;
elseif strcmp(sect,'preturn')
    int = 2;
elseif strcmp(sect,'turn')
    int = 3:4;
elseif strcmp(sect,'lat_arm')
    int = 5:6;
elseif strcmp(sect,'reward')
    int = 7:8;
elseif strcmp(sect,'delay')
    int = 9:12;
elseif strcmp(sect,'wheel')
    int = 13;
elseif strcmp(sect,'all')
    int = 1:13;
else
    error('Section is not a recognized. Options are [mid_arm, preturn, turn, lat_arm, reward, delay, wheel, all]');
end

if nXtra>0
    sect_times = varargin{1};
    sect_times = sect_times(int);
    varargout{1} = sect_times(sect_times>0);
end