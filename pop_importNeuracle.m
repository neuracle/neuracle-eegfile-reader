function [EEG, command] = pop_importNeuracle(filename,foldname, varargin);

% pop_importNeuracle() - import data files produced by Neuracle EEG Recorder into EEGLAB
% Usage:
%   >> [EEG, command] = pop_importNeuracle(); % pop up window
%   >> [EEG, command] = pop_importNeuracle( 'data.bdf', foldname );
%   >> [EEG, command] = pop_importNeuracle( {'data.bdf','evt.bdf'}, foldname, options );
%   >> [EEG, command] = pop_importNeuracle( {'data.bdf','data.1.bdf','data.2.bdf'}, foldname, options );
%
% Inputs:
%   filename - [string] one *.bdf file,
%               [cell]  more than one *.bdf files;
%   foldname - [string] folder name;

% Optional inputs:
%   'channels'   - [integer array] list of channel indices
%   'blockrange' - [min max] time limits for importing data
%   'memorymapped' - ['on'|'off'] import memory mapped file (useful if 
%                  encountering memory errors). Default is 'off'.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: read_bdf.m
%
% Author: Xiaoshan HUANG, hxs@neuracle.cn
%         Junying FANG, fangjunying@neuracle.cn
%
% Versions:
%    v1.0.0: 2017-09-27, orignal
%    v1.0.1: 2018-08-12, fix bug for shifting event latency for multiple data files
%    v1.0.2: 2019-03-28, correct 'boundary' of event and fread UTF-8 annotations
%    v1.1.0: 2019-12-26, add user panel used for selecting channels and blockrange of data
%    v1.1.1: 2020-11-13, add function input argument, 'filename'
% Copyright (c) 2020 Neuracle, Inc. All Rights Reserved. http://neuracle.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = [];
command = '';

if ~plugin_askinstall('NeuracleEEGFileReader', 'read_bdf'), return; end
n = 0;m = 0;
flagHEEG = 0;
if nargin < 1
	% ask user
    % --------
    ButtonName = questdlg2('Do you want to import a file or a folder?', ...
                           'Neuracle import', ...
                           'Folder', 'File', 'File');
    % collect data files and an evt file
    % ----------------------------------
    if strcmpi(ButtonName, 'file')
        [filename, foldname] = uigetfile({'*.bdf;*.edf';'*.*'}, 'Choose a file or header file -- pop_importNeuracle()','MultiSelect', 'on'); 
        drawnow;
        if ~iscell(filename),filename = {filename};end
        try if filename(1) == 0 ,return; end,  end
        for i = 1:length(filename)
            str = strsplit(filename{1,i},'.');
            if strcmpi(str{1,1},'data'), n= n+1;datafilename{1,n} = fullfile(foldname,filename{1,i});end
            if strcmpi(str{1,1},'evt'),m = m+1;evtfilename{1,m} = fullfile(foldname,filename{1,i});end
        end
    else
        foldname = uigetdir('*.*', 'Choose a folder -- pop_importNeuracle()'); 
        drawnow;
        if foldname(1) == 0 ,return; end
        sepfile = dir(foldname);names = {sepfile.name}';sepfile(ismember(names,{'.','..'})) = [];
        for i = 1:length(sepfile)
            if sepfile(i).isdir
                flagHEEG = 1;
                temp = dir(fullfile(foldname,sepfile(i).name));tempnames = {temp.name}';temp(ismember(tempnames,{'.','..'})) = [];
                for j = 1:length(temp)
                    str = strsplit(temp(j).name,'.');
                    if strcmpi(str{1,1},'data'),n= n+1; datafilename{1,n} = fullfile(foldname,sepfile(i).name,temp(j).name);end
                    if strcmpi(str{1,1},'evt'),m = m+1;evtfilename{1,m} = fullfile(foldname,sepfile(i).name,temp(j).name);end
                end
            else
                str = strsplit(sepfile(i).name,'.');
                if strcmpi(str{1,1},'evt')&strcmpi(str{end},'bdf'),m = m+1;evtfilename{1,m} = fullfile(foldname,sepfile(i).name);end
                if strcmpi(str{1,1},'data')& strcmpi(str{end},'bdf'),n = n+1;datafilename{1,n} = fullfile(foldname,sepfile(i).name);end   
            end
        end
    end
elseif  nargin == 1
    disp('wrong input');
else
    if isempty(foldname) || ~ischar(foldname)
         disp('wrong input');
    end     
    if isempty(filename)
        sepfile = dir(foldname);names = {sepfile.name}';sepfile(ismember(names,{'.','..'})) = [];
        for i = 1:length(sepfile)
            if sepfile(i).isdir
                flagHEEG = 1;
                temp = dir(fullfile(foldname,sepfile(i).name));tempnames = {temp.name}';temp(ismember(tempnames,{'.','..'})) = [];
                for j = 1:length(temp)
                    str = strsplit(temp(j).name,'.');
                    if strcmpi(str{1,1},'data'),n= n+1; datafilename{1,n} = fullfile(foldname,sepfile(i).name,temp(j).name);end
                    if strcmpi(str{1,1},'evt'),m = m+1;evtfilename{1,m} = fullfile(foldname,sepfile(i).name,temp(j).name);end
                end
            else
                str = strsplit(sepfile(i).name,'.');
                if strcmpi(str{1,1},'evt')&strcmpi(str{end},'bdf'),m = m+1;evtfilename{1,m} = fullfile(foldname,sepfile(i).name);end
                if strcmpi(str{1,1},'data')& strcmpi(str{end},'bdf'),n = n+1;datafilename{1,n} = fullfile(foldname,sepfile(i).name);end   
            end
        end
    else
        if ~iscell(filename),filename = {filename};end
        try if filename(1) == 0 ,return; end,  end
        for i = 1:length(filename)
            str = strsplit(filename{1,i},'.');
            if strcmpi(str{1,1},'data'), n= n+1;datafilename{1,n} = fullfile(foldname,filename{1,i});end
            if strcmpi(str{1,1},'evt'),m = m+1;evtfilename{1,m} = fullfile(foldname,filename{1,i});end
        end
    end 
end
    
if n == 0 
    disp('Current folder/files dose not contain data.bdf !');
    return;
end

if m > 1 
    disp('Only one evt.bdf file is support !');
    return;
end

% ---------------
% read file head to get infos
% ---------------------------
disp('Reading data file header...');
datafilelength = length(datafilename);
nTrials = 0;
datapnts = zeros(1, datafilelength);
T0 = [];
srate = [];
HDR =cell(1,datafilelength);
for i = 1:datafilelength
    hdr = read_bdf(datafilename{1,i});
    nTrials = nTrials + hdr.nTrials;
    datapnts(1,i) = hdr.nSamples;
    T0 = [T0; hdr.T0];
    srate = [srate,hdr.Fs];
    HDR{1,i} = hdr;
end
if length(unique(srate)) > 1, return,end
srate = unique(srate);
for j=1:size(T0,1)
    startsecs(j) = T0(j,3)*946080000+T0(j,2)*2592000+T0(j,1)*86400+T0(j,4)*3600+T0(j,5)*60+T0(j,6);
end
[startsecs,indexsorted] = sort(startsecs);
T0 = T0(indexsorted,:);
datafilename = datafilename(indexsorted); 
HDR = HDR(indexsorted);
if nargin < 1
    eeglab_options;
    uilist = { { 'style' 'text' 'String' 'Channel list (defaut all):' } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' [ 'Data range (in seconds) to read (default all [0 ' num2str(sum(datapnts)/srate) '])' ] } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'String' 'Extract event' } ...
                 { 'style' 'checkbox' 'string' '' 'value' m 'enable' 'on' } {} ...
                 { 'style' 'checkbox' 'String' 'Import as memory mapped file (use if out of memory error)' 'value' option_memmapdata } };
    geom = { [3 1] [3 1] [3 0.35 0.5] [1] };
    result = inputgui( geom, uilist, 'pophelp(''pop_importNeuracle'')', 'Load data using Neuracle-IO -- pop_importNeuracle()');
    if length(result) == 0 return; end
    % decode GUI params
    % -----------------

    options = {};
    try  
        eval( [ '[' result{1} ']' ]);
    catch
        chanlabels = lower(hdr.label);
        stridx = lower(strsplit(result{1},' '));
        [Lia,Locb] = ismember(stridx,chanlabels);
        result{1} = int2str(Locb(Lia)');
        if ~isempty(find(Lia == 0)),disp(['can not find channel: ' stridx{~Lia}]);end
    end
    if ~isempty(result{1}), options = { options{:} 'channels'   eval( [ '[' result{1} ']' ] ) }; end
    if ~isempty(result{2}), options = { options{:} 'blockrange' eval( [ '[' result{2} ']' ] ) }; end
    if length(result) > 2
        if ~result{3}, options = { options{:} 'importevent'  'off'  }; end
        if  result{4}, options = { options{:} 'memorymapped' 'on' }; end
    end
else
    options = varargin; 
end
% decode imput parameters
% -----------------------
g = struct(options{:});
if ~isfield(g, 'blockrange'), g.blockrange = []; end
if ~isfield(g, 'channels'), g.channels = []; end
if ~isfield(g, 'importevent'), if m == 0, g.importevent = 'off'; else, g.importevent = 'on'; end,end
if ~isfield(g, 'memorymapped'), g.memorymapped = 'off'; end

% import event
% -----------
comments = datafilename;
if strcmpi(g.importevent,'on')
    comments{end+1} = evtfilename{1,1};
    hdrEvt = read_bdf(evtfilename{1,1});
    evt = hdrEvt.event;
    evt = cell2mat(evt);
    event = struct('type',[],'latency',[]);
    if isstruct(evt)
        for i = 1:size(evt,2)
            event(i).type = evt(i).eventvalue;
            event(i).latency = round(evt(i).offset_in_sec*srate);
        end
    end
    Event = event;
    
    if datafilelength >1
        startpnts = floor((startsecs-startsecs(1))*srate);  %relative time
        endpnts = zeros(1,datafilelength);
        pausepnts = 0;
        Event = struct('type',[],'latency',[]);
        N = 0;
        datapnts = datapnts(indexsorted);
        for j = 1:datafilelength
            endpnts(j) = startpnts(j)+datapnts(j);
            if j > 1
                pausepnts = pausepnts+(startpnts(j)-endpnts(j-1));
            end
            if ~isempty([event(:).latency])
                for k = 1:length(event)
                    if startpnts(j) <= event(k).latency && endpnts(j) >= event(k).latency
                        event(k).latency = event(k).latency-pausepnts;
                        N = N+1;
                        Event(N) = event(k);
                    end
                end
            end
            if j < datafilelength
                eventBoundary = struct('type','boundary','latency',endpnts(j)-pausepnts+1); % relative ending
                N = N+1;
                Event(N) = eventBoundary;
            end

         end
        boundaries = findboundaries(Event);
        if boundaries(end)> sum(datapnts)
            Event(end) = [];
        end
    end
    
else
    if datafilelength > 1
        startpnts = floor((startsecs-startsecs(1))*srate);  %relative time
        endpnts = zeros(1,datafilelength);
        pausepnts = 0;
        Event = struct('type',[],'latency',[]);
        N = 0;
        datapnts = datapnts(indexsorted);
        for j = 1:datafilelength
            endpnts(j) = startpnts(j)+datapnts(j);
            if j > 1
                pausepnts = pausepnts+(startpnts(j)-endpnts(j-1));
            end
            if j < datafilelength
                eventBoundary = struct('type','boundary','latency',endpnts(j)-pausepnts+1); % relative ending
                N = N+1;
                Event(N) = eventBoundary;
            end
         end
        boundaries = findboundaries(Event);
        if boundaries(end)> sum(datapnts)
            Event(end) = [];
        end        
    else
        Event = struct('type',[],'latency',[]);
    end
end
% import data
% -----------
EEG = eeg_emptyset;
EEG.srate = srate;
EEG.trials = nTrials;
EEG.ref = 'common';
chaninfo.plotrad = [];
chaninfo.shrink = [];
chaninfo.nosedir = '+X';
chaninfo.icachansind = [];
EEG.chaninfo = chaninfo;
patientInfo = strsplit(hdr.orig.PID);
patientName = strtrim(patientInfo{4});
if ~strcmpi(patientName,'X')
    EEG.setname = [EEG.setname ' ' patientName ];     
end

% organize EEG info
% -------------
if length(g.blockrange) == 1,fprintf('Data range should contain two elements \n ');return,end
if isempty(g.blockrange)
    begsample = 1;
    endsample = sum(datapnts);
else
    begsample = floor(g.blockrange(1)*EEG.srate) + 1;
    endsample = floor(g.blockrange(end)*EEG.srate)-1;
end
if isempty(g.channels), chanidx = 1:hdr.nChans; else,chanidx = g.channels;end
EEG.chanlocs =  hdr.chanlocs(chanidx);
EEG.nbchan = length(chanidx);
EEG.pnts = endsample-begsample+1;
EEG.xmax = EEG.pnts/EEG.srate;
EEG.times = [1:EEG.pnts]*1000/EEG.srate; 

% read data
% ----------
fprintf('Reading data ...\n');

%%%%% only one data.bdf
if isempty([Event(:).latency]) & ~flagHEEG
    EEG.comments = char(comments);
    dat = read_bdf(datafilename{1,1}, HDR{1,1}, begsample, endsample, chanidx);
    if strcmpi(g.memorymapped,'on')
        alldata = mmo([], [EEG.nbchan,EEG.pnts]);
        for i = 1:EEG.nbchan
            alldata(i,:) = dat(i,:);
        end
        EEG.data = alldata;
    else
        EEG.data = dat;
    end
    EEG = eeg_checkset(EEG);
    % history
    % -------
    if ischar(foldname)
        if isempty(options)
            command = sprintf('EEG = pop_importNeuracle(''%s'');',foldname);
        else
            command = sprintf('EEG = pop_importNeuracle(''%s'', %s);', foldname, vararg2str(options));
        end
    end
    return
end

%%%%% more than one data.bdf files ,  HEEG
[~,I] = sort([Event(:).latency]);
Event = Event(I);
latency = [Event(:).latency];
[Lia,~] = ismember( {Event(:).type},'boundary');
boundarylatency = latency(Lia);
idxEvent =  find(begsample<= latency & latency <= endsample);
if ~isempty(idxEvent)
    for k = 1:length(idxEvent)
        EEG.event(k).type = Event(idxEvent(k)).type;
        EEG.event(k).latency = Event(idxEvent(k)).latency - begsample +1;
    end
end
 
if isempty(boundarylatency)
    dat = read_bdf(datafilename{1,1}, HDR{1,1}, begsample, endsample, chanidx);
    idx = 1;
else
    idxboundary = find (boundarylatency < endsample &  boundarylatency > begsample);
    if isempty(idxboundary) 
        idx = find(endsample < boundarylatency );
        if isempty(idx),idx = length(boundarylatency)+1;else,idx = idx(1);end
        if idx > 1,  begsample = begsample - boundarylatency(idx-1); endsample = endsample - boundarylatency(idx-1);end
        dat = read_bdf(datafilename{1,idx},HDR{1,idx},begsample,endsample,chanidx);
    else
        idx = unique([idxboundary,idxboundary+1]);
        for k = 1:length(idx)
            disp(['(' num2str(k) '/' num2str(length(idx)) '):' datafilename{1,idx(k)}]);
            if k == 1
                if idx(1) >1, newbegsample = begsample - boundarylatency(idx(1)-1); else;newbegsample = begsample;end
                newendsample = datapnts(idx(k));
            elseif k == length(idx)
                newbegsample = 1; newendsample = endsample-boundarylatency(idx(k)-1)+1;
            else
                newbegsample = 1;newendsample = datapnts(idx(k));
            end
            Data{k} = read_bdf(datafilename{1,idx(k)},HDR{1,idx(k)},newbegsample,newendsample,chanidx);
        end
        dat = cat(2,Data{:});
    end  
end 

if strcmpi(g.memorymapped,'on')
    alldata = mmo([], [EEG.nbchan,EEG.pnts]);
    for i = 1:EEG.nbchan
        alldata(i,:) = dat(i,:);
    end
    EEG.data = alldata;
else
    EEG.data = dat;
end
EEG.comments = char(comments{1,[idx,end]});
EEG = eeg_checkset(EEG);
% history
% -------
if ischar(foldname)
    if isempty(options)
        command = sprintf('EEG = pop_importNeuracle(''%s'');',foldname);
    else
        command = sprintf('EEG = pop_importNeuracle(''%s'', %s);', foldname, vararg2str(options));
    end
end
return
end

function boundaries = findboundaries(event)
if isfield(event, 'type') & isfield(event, 'latency') & cellfun('isclass', {event.type}, 'char')
    % Boundary event indices
    boundaries = strmatch('boundary', {event.type});
    % Boundary event latencies
    boundaries = [event(boundaries).latency];
    % Shift boundary events to epoch onset
    boundaries = fix(boundaries + 0.5);
    % Remove duplicate boundary events
    boundaries = unique(boundaries);
    % Epoch onset at first sample?
    if isempty(boundaries) || boundaries(1) ~= 1
        boundaries = [1 boundaries];
    end
else
    boundaries = 1;
end
end

