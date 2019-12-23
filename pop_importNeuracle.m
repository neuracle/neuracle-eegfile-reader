function EEG = pop_importNeuracle(pathname)
%Import EEG data files produced by Neuracle EEG Recorder into EEGLAB
% Syntax:  
%     EEG = readbdfdata(); % pop up window
%     EEG = readbdfdata(pathname) %read the data of all the channels
% 
% Inputs:
%     pathname: data files path
%
% Outputs:
%     EEG data structure
%
% Example:
%     EEG = readbdfdata(pathname)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: readbdfdata.m
%
% Author: Xiaoshan Huang, hxs@neuracle.cn
%         Junying FANG, fangjunying@neuracle.cn
%
% Versions:
%    v1.0: 2017-09-27, orignal
%    v1.1: 2018-08-12, update readbdfdata.m
% Copyright (c) 2017 Neuracle, Inc. All Rights Reserved. http://neuracle.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    [filename, pathname] = uigetfile({'*.bdf;*.edf';'*.*'}, 'Pick a recorded EEG data file','MultiSelect', 'on');
end

EEG = readbdfdata(filename, pathname);
eeglab redraw;
end

