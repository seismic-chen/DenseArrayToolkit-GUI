function varargout = seis_Find(varargin)
% SEIS_FIND: a M-file to fetch seismograms from IRIS. It is a GUI 
% program which requires seis_Find.fig
% 
% [event, seiswav, seisatt, channelinfo] = seis_Find(myEvent,mySeis)
% 
% Inputs:
%       myEvent:   an object of class MatEvent
%       myNetwork: an object of class MatNetwork
%       mySeis:    an object of class MatSeismogram
% Outputs:
%       event:       a structure array containing the tracking event 
%                    information
%       seiswav:     a cell array containing retrieved seismograms
%       seisatt:     a structure array saved attributions of retrieved
%                    seismograms
%       channelinfo: a structure array kept channel information
%
% Attention: 
% For multiple retrieves, the function return just reflects the
% last part. All retrieved seismograms and channel information
% are saved into a set of files named by user.
% 
% See also: MatEvent, MatSeismogram, get
 
% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Sep, 2003
%
% Editied by Darren Hart , on 9-25-2004
% Sandia National Laboratory, dhart@sanida.gov
% Added use of MatNetwork class. Origin ID instead of origin index.
% on 10-15-2004 for IRIS instrument response download beta, on 02-01-2004
% for NCEDC wfm fix.

% Edit the above text to modify the response to help seis_Find

% Last Modified by GUIDE v2.5 26-Feb-2004 21:06:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @seis_Find_OpeningFcn, ...
                   'gui_OutputFcn',  @seis_Find_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before seis_Find is made visible.
function seis_Find_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seis_Find (see VARARGIN)

% Choose default command line output for seis_Find
handles.output = hObject;
handles.seisgrampick=[];        % set variable for seismogram picking
handles.selection=[];
handles.seis=[];
handles.event=[];
handles.pick=[];
handles.save.seiswav=[];
handles.save.seisatt=[];
handles.save.channelinfo=[];
handles.save.network = [];
% Update handles structure
guidata(hObject, handles);

if nargin <5
    beep;
    msgbox('MatSeismogram and MatEvent objects required...','Remind','warn');
    return;
end

handles.myevent = varargin{1};
handles.mynetwork = varargin{2};
handles.myseis = varargin{3};
% Update handles structure
guidata(hObject, handles);

if ~isa(handles.myseis, 'MatSeismogram')
    beep;
    msgbox('MatSeismogram object required...','Remind','warn');
    return;
end

if ~isa(handles.myevent, 'MatEvent');
    beep;
    msgbox('MatEvent object required...','Remind','warn');
    return;
end
    
% disable other buttons
set(handles.retrieve_button,'Enable','off');
set(handles.select_button,'Enable','off');
set(handles.siteatt_button,'Enable','off');
set(handles.take_event,'Enable','off');
set(handles.save_button,'Enable','off');
set(handles.list_button,'Enable','off');
set(handles.map_button,'Enable','off');

cook=which('cookie.mat');
if ~exist(cook)
    set(handles.plot_button,'Enable','off');
    set(handles.total_event, 'String', 'No Events');
else
    seis=[];
    selection=[];
    events=[];
    % initialize
    
    load cookie;
    if exist('selection') && ~isempty(selection)
        handles.selection = selection;
    end
    if exist('seis') && isempty(seis) && isa(seis, 'MatSeismogram')
            handles.seis = seis;
    end
    if exist('events') && ~isempty(events) && isa(events, 'MatEvent')
        handles.event = events;
    end
    % Update handles structure
    guidata(hObject, handles);
    
    if exist('selection') && ~isempty(selection)
    % fill parameters for event
        set(handles.start_time,'String',handles.selection.origin.time);
        set(handles.event_long,'String',num2str(handles.selection.origin.longitude));
        set(handles.event_lat,'String',num2str(handles.selection.origin.latitude));
        set(handles.event_depth,'String',num2str(handles.selection.origin.depth));
    
    % fill parameters for site locating
        set(handles.network_edit, 'String', selection.p.network);
        tmp='';
        for i=1:length(selection.p.station)
            tmp=[tmp,' ',selection.p.station{i}];
        end
        set(handles.station_edit,'String',tmp);
        clear tmp;
        set(handles.site_edit, 'String', selection.p.site);
        set(handles.channel_edit, 'String', selection.p.channel);
    
    % fill parameters for selection
        switch lower(handles.selection.area.type)
            case 'boxarea'
                set(handles.select_area, 'Value', get(handles.select_area,'Max'));
                set(handles.select_dist, 'Value', get(handles.select_dist,'Min'));
                set(handles.area_longstart,'String',num2str(handles.selection.area.minLon));
                set(handles.area_longend,'String',num2str(handles.selection.area.maxLon));
                set(handles.area_latstart,'String',num2str(handles.selection.area.minLat));
                set(handles.area_latend,'String',num2str(handles.selection.area.maxLat));
            case 'distazim'
                set(handles.select_dist, 'Value', get(handles.select_dist,'Max'));
                set(handles.select_area, 'Value', get(handles.select_area,'Min'));
                set(handles.from_dist,'String',num2str(handles.selection.area.minDist));
                set(handles.to_dist,'String',num2str(handles.selection.area.maxDist));
                set(handles.from_azimuth,'String',num2str(handles.selection.area.minAzim));
                set(handles.to_azimuth,'String',num2str(handles.selection.area.maxAzim));
        end
    
    % fill phase parameters
        if ~isempty(handles.selection.phase)
            set(handles.phase_select, 'Value', get(handles.phase_select,'Max'));
            tmp=get(handles.phase1,'String');
            if strcmp('p,P,Pdiff,PKIKP', handles.selection.phase.fromPhase)
                handles.selection.phase.fromPhase='First P';
            end
            if strcmp('s,S,SKS,SKIKS', handles.selection.phase.fromPhase)
                handles.selection.phase.fromPhase='First S';
            end
            tmp=strcmp(tmp,handles.selection.phase.fromPhase);
            tmp=find(tmp==1);   % get the index of phase1
            set(handles.phase1,'Value',tmp);

            if strcmp('p,P,Pdiff,PKIKP', handles.selection.phase.toPhase)
                handles.selection.phase.toPhase='First P';
            end
            if strcmp('s,S,SKS,SKIKS', handles.selection.phase.toPhase)
                handles.selection.phase.toPhase='First S';
            end
            tmp=get(handles.phase2,'String');
            tmp=strcmp(tmp,handles.selection.phase.toPhase);
            tmp=find(tmp==1);   % get the index of phase2
            set(handles.phase2,'Value',tmp);
            clear tmp;
        
            tmp=ceil(handles.selection.phase.fromRelTime/60);   % set the relative time to phase1
            set(handles.from_time, 'String', num2str(tmp));
    
            tmp=ceil(handles.selection.phase.toRelTime/60);   % set the relative time to phase2
            set(handles.to_time, 'String', num2str(tmp));
            clear tmp;
        end
    
    end

    if exist('pick') && ~isempty(pick)
        handles.pick=pick;
        if length(pick) >1
            set(handles.take_event,'Enable','on');
            set(handles.take_event,'String',pick);
            set(handles.take_event,'Value',1.0);
        end
        set(handles.total_event, 'String',[num2str(length(pick)), ' Events']);
    else
        if exist('events') && ~isempty(events) && isa(events, 'MatEvent')
            set(handles.total_event, 'String',[num2str(handles.event.eventNum), ' Events']);
            if handles.event.eventNum >1
                set(handles.take_event,'Enable','on');
            end
            for kk=1:handles.event.eventNum;handles.pick{kk}=num2str(kk);end;
            set(handles.take_event,'String',handles.pick);
            set(handles.take_event,'Value',1.0);
            starting=handles.event.preferredOrigin(1);
            set(handles.start_time,'String',starting.originTime);
            set(handles.event_long,'String',starting.location.longitude);
            set(handles.event_lat,'String',starting.location.latitude);
            set(handles.event_depth,'String',starting.location.depth);
        end
    end
    
    if exist('seis') && ~isempty(seis) && isa(seis, 'MatSeismogram')
        set(handles.seis_no_end, 'String', num2str(get(handles.seis,'nseismogram')));
        set(handles.list_button,'Enable','on');
        set(handles.map_button,'Enable','on');
    end
    clear seis selection events;

    % Update handles structure
    guidata(hObject, handles);
end         
    
% initial load file functions
set(handles.loadFilename_text,'Visible','off');
set(handles.loadFilename_edit,'Visible','off');


% UIWAIT makes seis_Find wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = seis_Find_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles.selection)
    varargout{1} = handles.selection.origin;
else
    varargout{1} = handles.selection;
end

varargout{2} = handles.save.seiswav;
varargout{3} = handles.save.seisatt;
varargout{4} = handles.save.channelinfo;
varargout{5} = handles.save.network;

delete(handles.figure1);

% --- Executes during object creation, after setting all properties.
function start_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function start_time_Callback(hObject, eventdata, handles)
% hObject    handle to start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_time as text
%        str2double(get(hObject,'String')) returns contents of start_time as a double


% --- Executes on button press in EventFinder_button.
function EventFinder_button_Callback(hObject, eventdata, handles)
% hObject    handle to EventFinder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disable buttons
set(handles.retrieve_button,'Enable','off');
set(handles.select_button,'Enable','off');
set(handles.siteatt_button,'Enable','off');
set(handles.plot_button,'Enable','off');
set(handles.list_button,'Enable','off');
set(handles.map_button,'Enable','off');

handles.pick=[];
% Update handles structure
guidata(hObject, handles);
% clear events

[handles.param, handles.event, pick]=event_Find(handles.myevent);
preforids = get(handles.event,'preferredOrigin');

if isempty(pick)
    beep;
    set(handles.message,'String','Select event first...');
else
    for k=1:length(pick)
        %handles.pick{k}= num2str(pick(k));
        handles.pick{k}= num2str(preforids(k).id);
    end;
    % Update handles structure
    guidata(hObject, handles);

    starting=handles.event.preferredOrigin(pick(1));
    set(handles.start_time,'String',starting.originTime);
    set(handles.event_long,'String',starting.location.longitude);
    set(handles.event_lat,'String',starting.location.latitude);
    set(handles.event_depth,'String',starting.location.depth);
    if length(pick) >1
        set(handles.take_event,'Enable','on');
    end
    set(handles.take_event,'String',handles.pick);
    set(handles.take_event,'Value',1.0);
end
set(handles.total_event, 'String', [num2str(length(handles.pick)),' Events']);
set(handles.message, 'String', [num2str(length(handles.pick)),' Events']);


% --- Executes on button press in check_button.
function check_button_Callback(hObject, eventdata, handles)
% hObject    handle to check_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% display message
set(handles.message, 'String', 'Searching available seismograms...');

% disable other buttons
set(handles.retrieve_button,'Enable','off');
set(handles.select_button,'Enable','off');
set(handles.siteatt_button,'Enable','off');
set(handles.plot_button,'Enable','off');
set(handles.EventFinder_button,'Enable','off');
set(handles.list_button,'Enable','off');
set(handles.map_button,'Enable','off');


% get all parameters required for checking seismograms
st=get(handles.start_time,'String');
duration_start=str2num(get(handles.duration_start,'String'))*60;

% get the starting time
ending=timeshift(st, duration_start);
handles.p.time.startTime=ending;

% get the ending time
duration_end=str2num(get(handles.duration_end,'String'))*60;
ending=timeshift(st, duration_end);
handles.p.time.endTime=ending;
% set time window

handles.p.network=upper(get(handles.network_edit,'String'));
% get network ID

station=get(handles.station_edit,'String');
i=1;
while (any(station))
  [chopped,station] = strtok(station);
% here, delimiter is whitespace
  handles.p.station{i} = upper(chopped);
  i=i+1;
end
% get station ID

handles.p.site=upper(get(handles.site_edit,'String'));
% get site ID

handles.p.channel=upper(get(handles.channel_edit,'String'));

% Update handles structure
guidata(hObject, handles);

% query about available seismograms
[handles.seis, nseis]=seis_query(handles.myseis,handles.p.time,handles.p.network,...
    handles.p.station,handles.p.site,handles.p.channel);
handles.nseis=nseis;
handles.gui.seis = handles.seis;
if nseis<1
    beep;
    set(handles.select_button,'Enable','off');
    set(handles.message,'String', 'No seismograms found, run query again');
else
    % find sites with seismograms
    handles.seis = channel_location(handles.seis);
    % Update handles structure
    guidata(hObject, handles);
    
    % display done message
    set(handles.message,'String',[num2str(nseis), ' seismograms found']);
    set(handles.select_button,'Enable','on');
end

% activate buttons
set(handles.EventFinder_button,'Enable','on');
%handles.mynetwork = retrieve_network(handles.mynetwork);
%handles.mynetwork = choose_network(handles.mynetwork,'code',handles.p.network);
%chosenNet = handles.mynetwork.chosenNetwork;
%handles.mynetwork = set_current_network(handles.mynetwork,chosenNet(1));
%handles.mynetwork = find_station(handles.mynetwork,'area', [-90 90 -180 180] );
%handles.mynetwork = find_station_channel(handles.mynetwork, handles.p.channel,[]);
%handles.mynetwork = find_channel_by_code(handles.mynetwork, handles.p.station,handles.p.channel,[]);
%handles.mynetwork = station_pool_add(handles.mynetwork);
%[handles.mynetwork stachanindex] = get_response(handles.mynetwork, 'all', 'all',time2string(handles.p.time.startTime));

%count = 0;
%for i=1:length(stachanindex.staIndex)
%    network.station = handles.mynetwork.networkstation(stachanindex.staIndex(i));
%    for j=1:length(stachanindex.chanIndex)
%        count = count +1;
%        handles.network.instrument(count) = network.station.channel(stachanindex(i).chanIndex(j)).instrumentation;
%    end
%end
% Update handles structure
%guidata(hObject, handles);
%display('instrument response files downloaded...')
% --- Executes on button press in select_button.
function select_button_Callback(hObject, eventdata, handles)
% hObject    handle to select_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.message, 'String', 'Selecting seismograms...');

% set event info
origin.time = get(handles.start_time,'String');
origin.longitude = str2num(get(handles.event_long,'String'));
origin.latitude = str2num(get(handles.event_lat,'String'));
origin.depth = str2num(get(handles.event_depth,'String'));
o_value = get(handles.take_event,'Value');
orid = get(handles.take_event,'String');
if ischar(orid)
    origin.orid = str2num(orid(o_value,:));
elseif iscell(orid)
    origin.orid = str2num(orid{o_value});
end
origin.orid = str2num(orid{o_value});
orig = get(handles.event,'preferredOrigin');
origin_temp = orig(o_value);
origin.mag_value = origin_temp.magnitude.value;
origin.mag_type = origin_temp.magnitude.type;
origin.catalog = origin_temp.catalog;
origin.contributor = origin_temp.contributor;

% get window information
win_start = timeshift( origin.time, str2num( get(handles.duration_start ,  'String') )*60 );
win_end =   timeshift( origin.time, str2num( get(handles.duration_end , 'String') )*60 );
time.fromRelTime = win_start;
time.toRelTime = win_end;
% select seismograms by a box area
if get(handles.select_area, 'Value') == get(handles.select_area, 'Max')
    area.type = 'boxArea';
    area.minLon=str2num(get(handles.area_longstart,'String'));
    area.maxLon=str2num(get(handles.area_longend,'String'));
    area.minLat=str2num(get(handles.area_latstart,'String'));
    area.maxLat=str2num(get(handles.area_latend,'String'));
end
    
% select seismograms by a azimuth distance
if get(handles.select_dist, 'Value') == get(handles.select_dist, 'Max')
    area.type = 'distAzim';
    area.minDist=str2num(get(handles.from_dist,'String'));
    area.maxDist=str2num(get(handles.to_dist,'String'));
    area.minAzim=str2num(get(handles.from_azimuth,'String'));
    area.maxAzim=str2num(get(handles.to_azimuth,'String'));
end

% select time window by phases
if get(handles.phase_select, 'Value') == get(handles.phase_select, 'Max')
    allphase1 = get(handles.phase1,'String');
    id1 = get(handles.phase1,'Value');
    phase.fromPhase=allphase1{id1};
    
    % define the first P and first S selection
    if id1 == 1
        phase.fromPhase = 'p,P,Pdiff,PKIKP';
    elseif id1==9
        phase.fromPhase = 's,S,SKS,SKIKS';
    end
    phase.fromRelTime = str2num(get(handles.from_time,'String'))*60;
    
    allphase2 = get(handles.phase2,'String');
    id2 = get(handles.phase2,'Value');
    phase.toPhase = allphase2{id2};
    
    % define the first P and first S selection
    if id2 == 1
        phase.toPhase = 'p,P,Pdiff,PKIKP';
    elseif id1==9
        phase.toPhase = 's,S,SKS,SKIKS';
    end

    phase.toRelTime = str2num(get(handles.to_time, 'String'))*60;
else
    phase=[];
    %phase.fromPhase = 'p,P,Pdiff,PKIKP';
    %phase.toPhase = 'p,P,Pdiff,PKIKP';
end

% select seismograms
[handles.seis, nselect] = seis_select(handles.seis, origin, area, phase, time, handles);


% backup parameters
handles.selection.p=handles.p;
handles.selection.origin=origin;
handles.selection.area=area;
handles.selection.phase=phase;
handles.nselect=nselect;

clear handles.mySeis;
clear handles.myEvent;

% Update handles structure
guidata(hObject, handles);

set(handles.retrieve_button,'Enable','on');
set(handles.return_button,'Enable','on');
set(handles.message,'String',[num2str(handles.nseis), ' seismograms found, ', num2str(handles.nselect), ' seismograms selected']);

% --- Executes on button press in retrieve_button.
function retrieve_button_Callback(hObject, eventdata, handles)
% hObject    handle to retrieve_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set event info

% display retrieve message
set(handles.message,'String','Retrieving seismograms...');

% get filename for saving data and the number of seismograms on each
% retrieving
[t_filen, nRetrieve]=getFilename(handles.selection.origin.time);

% keep filename and the number of retrievings so that when user presses
% save_button button to merge them together
handles.filename=strrep(t_filen,':','-');;
if nRetrieve~=0
    handles.seis=enable_iterator(handles.seis);
    handles.seis=init_iterator(handles.seis,nRetrieve);
    nIterat=ceil(handles.nselect/nRetrieve);
else
    nIterat=1;
    handles.seis=disable_iterator(handles.seis);
end

event=handles.selection.origin;
nTotal=0;   % Total number of retrieved seismograms
AllFiles=[];        % All temporary filenames
for i=1: nIterat
    t_filename=[t_filen,':',num2str(i)];
    t_filename=strrep(t_filename,':','-');
    if i==1
        set(handles.selectFiles,'String',t_filename);
    end
    AllFiles{i}=t_filename;
    % retrieve seismograms and channel info
    [handles.seis, nSeiswav]=seis_retrieve(handles.seis, 'selected');
    handles.seis=channel_info(handles.seis);
    % getting seismograms
    seiswav=get(handles.seis, 'seismogram');
    % getting seismogram attributes
    seisatt=get(handles.seis,'seismogramAttr');
    channelinfo = get(handles.seis,'seisChanInfo');
    save(t_filename, 'event','seiswav', 'seisatt', 'channelinfo');
    nTotal=nTotal + nSeiswav;
    set(handles.message,'String',[num2str(nTotal), ' Seismograms retrieved, ', num2str(i),...
        'th iteration: ', num2str(nSeiswav)]);
end

set(handles.selectFiles, 'String',AllFiles);    % show all temporary files
handles.save.seiswav=seiswav;
handles.save.seisatt=seisatt;
handles.save.channelinfo=channelinfo;
handles.save.network = handles.mynetwork;
handles.AllFiles=AllFiles;
% Update handles structure
guidata(hObject, handles);

% activate return button
set(handles.return_button,'Enable','on');
set(handles.plot_button,'Enable','off');
set(handles.save_button,'Enable','on');
%set(handles.siteatt_button,'Enable','on');
set(handles.list_button,'Enable','off');
set(handles.map_button,'Enable','off');

set(handles.seis_no_end,'String',num2str(nTotal));
isparamerr(handles.seis_no_start,1,num2str(nTotal),handles);
% display task finished message
set(handles.message,'String',[num2str(nTotal), ' Seismograms retrieved, ', num2str(nIterat),...
                'th iteration: ', num2str(nSeiswav),'. Retrieving Done']);

% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% display message
set(handles.message,'String','Ploting seismograms...');

% get the number of seismograms
nseis=length(handles.save.seiswav);
seiswav = handles.save.seiswav;
seisatt=handles.save.seisatt;
seischaninfo = handles.save.channelinfo;

pick=handles.seisgrampick;
% get seismogram selection

startn=str2num(get(handles.seis_no_start, 'String'));
endn=str2num(get(handles.seis_no_end, 'String'));
if startn > endn
     tt=endn;
     endn=startn;
     startn=tt;
end
% swap

if get(handles.one_by_one, 'Value') == get(handles.one_by_one, 'Max')
    % get settings for figures
    rows=str2num(get(handles.rows,'String'));
    cols=str2num(get(handles.columns,'String'));
    tols=rows*cols;

    % check if seis_find starting index is valid
    isparamerr(handles.seis_no_start,1,nseis,handles);

    % check if seis_find ending index is valid
    isparamerr(handles.seis_no_end,1,nseis,handles);

    % set the name of figure
    hseis=figure(1000);

    for ii=startn:endn
        if exist('hseis')
            figure(hseis);      % set the current figure
        else
            hseis=figure(1000);
        end
    
        if mod(ii-startn+1,tols)==1 
            clf;
            if ii+tols-1<endn
                endseisn=num2str(ii+tols-1);
            else
                endseisn=num2str(endn);
            end % the index number of the last seismogram in the figure
            set(hseis,'NumberTitle', 'off', 'Name', ['Plotting from SeismoGrams ', ...
                get(handles.seis_no_start, 'String'), ...
                ' to ', get(handles.seis_no_end, 'String'), ' currently from ', num2str(ii), ' to ',endseisn]);
            
            subplot(rows,cols,1)
        elseif mod(ii-startn+1,tols)==0
            subplot(rows,cols,tols);
        else
            subplot(rows,cols,mod(ii-startn+1,tols));
        end
        
        if ~isempty(pick)
            attr=seisatt(pick(ii));
            chan=seischaninfo(pick(ii));
            waveform=seiswav(pick(ii));
        else
            attr=seisatt(ii);
            chan=seischaninfo(ii);
            waveform=seiswav(ii);
        end
        
        waveform=waveform{1};
        sampleNum = attr.nPoint;
    
        % absolute time of the first sample
        seisStartTime = attr.beginTime;
    
        tx= timedif(handles.selection.origin.time,seisStartTime) + (0:sampleNum-1)*attr.sampleIntv*0.001;
        plot(tx,waveform);
        
        [delta,azim,backazim]=delaz(handles.selection.origin.latitude,handles.selection.origin.longitude, ...
            chan.location.latitude,  chan.location.longitude, 0);

        if mod(ii-startn+1,tols)==1 || tols==1
            title({[handles.selection.origin.time, '   Long:',num2str(handles.selection.origin.longitude),...
                        '\circ   Lat:',num2str(handles.selection.origin.latitude), '\circ   Dep:', ...
                        num2str(handles.selection.origin.depth),'km'];[attr.network, '.', ...
                        attr.station, '.', attr.site, '.',attr.channel,...
                        '   \Delta:', num2str(round(delta)),'\circ   Az:',num2str(round(azim)),'\circ']});
        else
            title([seisatt(ii).network, '.', ...
                        attr.station, '.', attr.site, '.',attr.channel,...
                        '   \Delta:', num2str(round(delta)),'\circ   Az:',num2str(round(azim)),'\circ']);
        end
            
        xlabel('Time after orginTime(s)');
    
    % pause
        if mod(ii-startn+1,tols)==0 & (ii-startn+1)~=(endn-startn+1)
            set(handles.message,'String','Press any key to continue plotting...');
            pause;
        end
    end;
    % display message
    set(handles.message,'String','Seismograms plotting done');
end

if get(handles.all_in_one, 'Value') == get(handles.all_in_one, 'Max')
    if exist('hseis')
        figure(hseis);      % set the current figure
    else
        hseis=figure(1000);
    end
    clf;
    set(hseis,'NumberTitle', 'off', 'Name', 'Plotting SeismoGrams ');
   
    hold off
    scal=str2num(get(handles.scale,'String'));

    nseis=abs(endn-startn+1);
    delta=zeros(nseis,1);azim=delta;backazim=delta;max_amp=delta;
    for k=startn:endn
        
        if isempty(pick)
            seis=seiswav{k};
            chan=seischaninfo(k);
            att=seisatt(k);
        else
            seis=seiswav{pick(k)};
            chan=seischaninfo(pick(k));
            att=seisatt(pick(k));
        end

        [delta(k),azim(k),backazim(k)]=delaz(...
                handles.selection.origin.latitude,handles.selection.origin.longitude, ...
                chan.location.latitude,  chan.location.longitude, 0);
        max_amp(k) = max(abs(seis));
        
        tdiff = timedif(handles.selection.origin.time,att.beginTime);
        tim = tdiff + [0:att.nPoint-1]*att.sampleIntv*.001;
        plot(tim, delta(k) + seis/max_amp(k)*scal)
        hold on
        ht=text(tdiff, delta(k),chan.stationCode);
        set(ht, 'Color',[1 0 0], 'FontAngle', 'italic','FontName','FixedWidth', 'FontSize',12, 'HorizontalAlignment', 'Right');
    end
    tp1=get(handles.plotphase,'Value');
    allphases = get(handles.plotphase,'String');
    tphase=allphases{tp1};
    
    % define the first P and first S selection
    if tp1 == 1
        tphase = 'p,P,Pdiff,PKIKP';
    elseif tp1==9
        tphase = 's,S,SKS,SKIKS';
    end
    
    tt=taupCurve('ak135',handles.selection.origin.depth,tphase);
	for k=1:length(tt); plot(tt(k).time,tt(k).distance,'-r');end

    %del=sort(delta)';TT=get_ttt(['P';'S'],handles.selection.origin.depth,del);plot(TT,del)
    ylabel('Distance (\circ)');
    xlabel('Time after origin time (s)')
    % display message
    set(handles.message,'String','Seismograms plotting done');
end

% --- Executes on button press in return_button.
function return_button_Callback(hObject, eventdata, handles)
% hObject    handle to return_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cook=which('cookie.mat');
if isempty(cook)
    cook='cookie.mat';
end
cook=cook(1:length(cook)-4);
pick=handles.pick;
selection=handles.selection;
events = handles.event;
save(cook, 'selection', 'events', 'pick');
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function network_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to network_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function network_edit_Callback(hObject, eventdata, handles)
% hObject    handle to network_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of network_edit as text
%        str2double(get(hObject,'String')) returns contents of network_edit as a double


% --- Executes during object creation, after setting all properties.
function station_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to station_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function station_edit_Callback(hObject, eventdata, handles)
% hObject    handle to station_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of station_edit as text
%        str2double(get(hObject,'String')) returns contents of station_edit as a double


% --- Executes during object creation, after setting all properties.
function site_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to site_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function site_edit_Callback(hObject, eventdata, handles)
% hObject    handle to site_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of site_edit as text
%        str2double(get(hObject,'String')) returns contents of site_edit as a double


% --- Executes during object creation, after setting all properties.
function channel_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function channel_edit_Callback(hObject, eventdata, handles)
% hObject    handle to channel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel_edit as text
%        str2double(get(hObject,'String')) returns contents of channel_edit as a double


% --- Executes on button press in siteatt_button.
function siteatt_button_Callback(hObject, eventdata, handles)
% hObject    handle to siteatt_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message,'String','Retrieving site information takes time, be patient...');
Total=0;
for i=1: length(handles.AllFiles);
    load(handles.AllFiles{i});
    for j=1: length(seisatt)
        satt(Total+j)=seisatt(j);
    end
    Total=Total+length(seisatt);
end
a0=datestr(now,30);
SiteattRetriever(satt);
a1=datestr(now,30);
a3=num2str(ceil(timedif(a0,a1)/60));
set(handles.message,'String',['Site-Attributes-Retrieving done in ',a3,' Mins']);

% --- Executes on button press in select_area.
function select_area_Callback(hObject, eventdata, handles)
% hObject    handle to select_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_area
if get(handles.select_area, 'Value') == get(handles.select_area, 'Max')
   set(handles.select_dist, 'Value', get(handles.select_dist,'Min'));
else
    set(handles.select_dist, 'Value', get(handles.select_dist,'Max'));
end

% --- Executes during object creation, after setting all properties.
function area_latstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to area_latstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function area_latstart_Callback(hObject, eventdata, handles)
% hObject    handle to area_latstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of area_latstart as text
%        str2double(get(hObject,'String')) returns contents of area_latstart as a double
isparamerr(hObject,-90,90,handles);

% --- Executes during object creation, after setting all properties.
function area_longend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to area_longend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function area_longend_Callback(hObject, eventdata, handles)
% hObject    handle to area_longend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of area_longend as text
%        str2double(get(hObject,'String')) returns contents of area_longend as a double
isparamerr(hObject, -180, 180, handles);

% --- Executes during object creation, after setting all properties.
function area_longstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to area_longstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function area_longstart_Callback(hObject, eventdata, handles)
% hObject    handle to area_longstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of area_longstart as text
%        str2double(get(hObject,'String')) returns contents of area_longstart as a double
isparamerr(hObject,-180,180,handles);

% --- Executes during object creation, after setting all properties.
function area_latend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to area_latend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function area_latend_Callback(hObject, eventdata, handles)
% hObject    handle to area_latend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of area_latend as text
%        str2double(get(hObject,'String')) returns contents of area_latend as a double
isparamerr(hObject,-90,90,handles);

% --- Executes on button press in select_dist.
function select_dist_Callback(hObject, eventdata, handles)
% hObject    handle to select_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_dist
if get(handles.select_dist,'Value') == get(handles.select_dist,'Max')
    set(handles.select_area, 'Value', get(handles.select_area, 'Min'));
else
    set(handles.select_area, 'Value', get(handles.select_area, 'Max'));
end

% --- Executes during object creation, after setting all properties.
function from_azimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function from_azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to from_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from_azimuth as text
%        str2double(get(hObject,'String')) returns contents of from_azimuth as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;


% --- Executes during object creation, after setting all properties.
function to_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function to_dist_Callback(hObject, eventdata, handles)
% hObject    handle to to_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to_dist as text
%        str2double(get(hObject,'String')) returns contents of to_dist as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;

% --- Executes during object creation, after setting all properties.
function from_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function from_dist_Callback(hObject, eventdata, handles)
% hObject    handle to from_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from_dist as text
%        str2double(get(hObject,'String')) returns contents of from_dist as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;

% --- Executes during object creation, after setting all properties.
function to_azimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function to_azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to to_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to_azimuth as text
%        str2double(get(hObject,'String')) returns contents of to_azimuth as a double
isparamerr(hObject, 0, 360, handles);

% --- Executes on button press in phase_select.
function phase_select_Callback(hObject, eventdata, handles)
% hObject    handle to phase_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of phase_select


% --- Executes during object creation, after setting all properties.
function duration_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to duration_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function duration_start_Callback(hObject, eventdata, handles)
% hObject    handle to duration_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of duration_start as text
%        str2double(get(hObject,'String')) returns contents of duration_start as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;


% --- Executes during object creation, after setting all properties.
function duration_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to duration_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function duration_end_Callback(hObject, eventdata, handles)
% hObject    handle to duration_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of duration_end as text
%        str2double(get(hObject,'String')) returns contents of duration_end as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;


% --- Executes during object creation, after setting all properties.
function from_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function from_time_Callback(hObject, eventdata, handles)
% hObject    handle to from_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from_time as text
%        str2double(get(hObject,'String')) returns contents of from_time as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;


% --- Executes during object creation, after setting all properties.
function to_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function to_time_Callback(hObject, eventdata, handles)
% hObject    handle to to_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to_time as text
%        str2double(get(hObject,'String')) returns contents of to_time as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;


% --- Executes during object creation, after setting all properties.
function phase1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in phase1.
function phase1_Callback(hObject, eventdata, handles)
% hObject    handle to phase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns phase1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from phase1


% --- Executes during object creation, after setting all properties.
function rows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function rows_Callback(hObject, eventdata, handles)
% hObject    handle to rows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rows as text
%        str2double(get(hObject,'String')) returns contents of rows as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');
end
clear tmp;

% --- Executes during object creation, after setting all properties.
function columns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to columns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function columns_Callback(hObject, eventdata, handles)
% hObject    handle to columns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of columns as text
%        str2double(get(hObject,'String')) returns contents of columns as a double
tmp=str2double(get(hObject,'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    set(handles.message, 'String','');
    set(hObject, 'ForegroundColor','black');  
end
clear tmp;

% --- Executes during object creation, after setting all properties.
function seis_no_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seis_no_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function seis_no_start_Callback(hObject, eventdata, handles)
% hObject    handle to seis_no_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seis_no_start as text
%        str2double(get(hObject,'String')) returns contents of seis_no_start as a double
pick=handles.seisgrampick;
if ~isempty(pick)
    ulim=length(pick);
else
    ulim=length(handles.seiswav);
end
isparamerr(hObject, 1,ulim ,handles);

% --- Executes during object creation, after setting all properties.
function seis_no_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seis_no_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function seis_no_end_Callback(hObject, eventdata, handles)
% hObject    handle to seis_no_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seis_no_end as text
%        str2double(get(hObject,'String')) returns contents of seis_no_end as a double
pick=handles.seisgrampick;
if ~isempty(pick)
    ulim=length(pick);
else
    ulim=length(handles.seiswav);
end
isparamerr(hObject, 1,ulim ,handles);
    

function isparamerr(hObject, lmit, umit, handles)
% Usage:    bObject,    the handle of one GUI edit component
%           umit,       the number's upper limit
%           lmit,       the number's lower limit
%           handles,    the handle struction

tmp=str2double(get(hObject, 'String'));
if isnan(tmp)
    beep; set(handles.message, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    if tmp > umit 
        beep;
        set(hObject, 'String', num2str(umit));
        set(hObject, 'ForegroundColor','red');
        set(handles.message, 'String', ['Must be less than ', num2str(umit), ' !']);
    elseif tmp <lmit
        beep;
        set(hObject, 'String', num2str(lmit));
        set(hObject, 'ForegroundColor','red');
        set(handles.message, 'String', ['Must be greater than ', num2str(lmit), ' !']);
    else
        set(hObject, 'ForegroundColor','black');
        set(handles.message, 'String', '');
    end
end
clear tmp;


% --- Executes during object creation, after setting all properties.
function event_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function event_name_Callback(hObject, eventdata, handles)
% hObject    handle to event_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_name as text
%        str2double(get(hObject,'String')) returns contents of event_name as a double


% --- Executes during object creation, after setting all properties.
function event_long_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_long (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function event_long_Callback(hObject, eventdata, handles)
% hObject    handle to event_long (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_long as text
%        str2double(get(hObject,'String')) returns contents of event_long as a double


% --- Executes during object creation, after setting all properties.
function event_lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function event_lat_Callback(hObject, eventdata, handles)
% hObject    handle to event_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_lat as text
%        str2double(get(hObject,'String')) returns contents of event_lat as a double


% --- Executes during object creation, after setting all properties.
function event_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function event_depth_Callback(hObject, eventdata, handles)
% hObject    handle to event_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_depth as text
%        str2double(get(hObject,'String')) returns contents of event_depth as a double


% --- Executes during object creation, after setting all properties.
function phase2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in phase2.
function phase2_Callback(hObject, eventdata, handles)
% hObject    handle to phase2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns phase2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from phase2


% --- Executes during object creation, after setting all properties.
function take_event_CreateFcn(hObject, eventdata, handles)
% hObject    handle to take_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function take_event_Callback(hObject, eventdata, handles)
% hObject    handle to take_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of take_event as text
%        str2double(get(hObject,'String')) returns contents of take_event as a double
idx=round(get(hObject,'Value'));
if length(handles.pick)>0
    %starting=handles.event.preferredOrigin(str2num(handles.pick{idx}));
    starting=handles.event.preferredOrigin(idx);
    set(handles.start_time,'String',starting.originTime);
    set(handles.event_long,'String',starting.location.longitude);
    set(handles.event_lat,'String',starting.location.latitude);
    set(handles.event_depth,'String',starting.location.depth);
else
    set(handles.message,'String', 'No event selected...');
end

% --- Executes on button press in all_in_one.
function all_in_one_Callback(hObject, eventdata, handles)
% hObject    handle to all_in_one (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of all_in_one

if get(handles.all_in_one, 'Value') == get(handles.all_in_one, 'Max')
    set(handles.one_by_one, 'Value', get(handles.one_by_one,'Min'));
else
    set(handles.one_by_one, 'Value', get(handles.one_by_one,'Max'));
end


% --- Executes during object creation, after setting all properties.
function scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function scale_Callback(hObject, eventdata, handles)
% hObject    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale as text
%        str2double(get(hObject,'String')) returns contents of scale as a double


% --- Executes on button press in one_by_one.
function one_by_one_Callback(hObject, eventdata, handles)
% hObject    handle to one_by_one (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of one_by_one

if get(handles.one_by_one, 'Value') == get(handles.one_by_one, 'Max')
    set(handles.all_in_one, 'Value', get(handles.all_in_one,'Min'));
else
    set(handles.all_in_one, 'Value', get(handles.all_in_one,'Max'));
end


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message, 'String', 'Saving all temporary files into one file, Please wait....');
clear handles.mySeis;
clear handles.myEvent;
clear handles.seis;

if find(strcmp(handles.AllFiles,handles.filename)==1)
    set(handles.message,'String',[handles.filename,' already exists, no need to merge']);
    return;
end

Total=0;
for i=1: length(handles.AllFiles);
    
    load(handles.AllFiles{i});
    for j=1: length(seisatt)
        swav(Total+j)=seiswav(j);
        satt(Total+j)=seisatt(j);
        scha(Total+j)=channelinfo(j);
    end
    Total=Total+length(seisatt);
end
event=handles.selection.origin;
seiswav=swav;
seisatt=satt;
channelinfo=scha;
save(handles.filename,'event','seiswav','seisatt','channelinfo'); 
set(handles.message, 'String', 'Saving done...');

% add the merged file to the popmenu for file-seclecting
allfiles=get(handles.selectFiles,'String');
allfiles{length(allfiles)+1}=handles.filename;
set(handles.selectFiles,'String',allfiles);
set(handles.selectFiles,'Value',length(allfiles));

clear seiswav;
clear seisatt;
clear channelinfo;

handles.AllFiles=allfiles;
handles.save.seiswav=swav;
handles.save.seisatt=satt;
handles.save.channelinfo=scha;
% Update handles structure
guidata(hObject, handles);

% activate buttons
set(handles.plot_button,'Enable','on');
set(handles.list_button,'Enable','on');
set(handles.map_button,'Enable','on');
%set(handles.save_button,'Enable','off');
set(handles.seis_no_end,'String',num2str(length(swav)));
isparamerr(handles.seis_no_start,1,num2str(length(swav)),handles);


% --- Executes during object creation, after setting all properties.
function plotphase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in plotphase.
function plotphase_Callback(hObject, eventdata, handles)
% hObject    handle to plotphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plotphase contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotphase


% --- Executes on button press in list_button.
function list_button_Callback(hObject, eventdata, handles)
% hObject    handle to list_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.message, 'String', 'Listing seismograms...');
handles.seisgrampick=[];    % initialization
handles.seisgrampick=SeisLister(handles.save, handles.selection.origin);
% Update handles structure
guidata(hObject, handles);

if ~isempty(handles.seisgrampick)
    set(handles.seis_no_end, 'String', num2str(length(handles.seisgrampick)));
else
    set(handles.seis_no_end, 'String', num2str(length(handles.save.seisatt)));
end
set(handles.message, 'String', 'Seismogram listing done');


% --- Executes on button press in map_button.
function map_button_Callback(hObject, eventdata, handles)
% hObject    handle to map_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message, 'String', 'Mapping seismogram stations...');
SeisMaper('Seisatt',handles.save.seisatt, ...
    'Channelinfo',handles.save.channelinfo,...
    'Event', handles.selection.origin, 'Selection',handles.seisgrampick);
set(handles.message, 'String', 'Station mapping done');


% --- Executes during object creation, after setting all properties.
function selectFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in selectFiles.
function selectFiles_Callback(hObject, eventdata, handles)
% hObject    handle to selectFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectFiles
set(handles.message,'String', 'Select one of the files for retrieved seismograms...');
Ind=get(handles.selectFiles, 'Value');
nameFile=handles.AllFiles{Ind};             % get the user-selected file
load(nameFile);
handles.save.seiswav=seiswav;
handles.save.seisatt=seisatt;
handles.save.channelinfo=channelinfo;
handles.seisgrampick=[1:length(seiswav)];
% Update handles structure
guidata(hObject, handles);

% activate buttons
set(handles.plot_button,'Enable','on');
set(handles.save_button,'Enable','on');
set(handles.list_button,'Enable','on');
set(handles.map_button,'Enable','on');
set(handles.seis_no_end,'String',num2str(length(seiswav)));
isparamerr(handles.seis_no_start,1,num2str(length(seiswav)),handles);


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.message,'String','Loading seimograms from a given file');

% activate load file functions
set(handles.loadFilename_text,'Visible','on');
set(handles.loadFilename_edit,'Visible','on');
%set(handles.siteatt_button,'Enable','on');
set(handles.save_button,'Enable','off');


% --- Executes during object creation, after setting all properties.
function loadFilename_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadFilename_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function loadFilename_edit_Callback(hObject, eventdata, handles)
% hObject    handle to loadFilename_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadFilename_edit as text
%        str2double(get(hObject,'String')) returns contents of loadFilename_edit as a double
% initial load file functions

filename=get(hObject,'String');
set(handles.loadFilename_text,'Visible','off');
set(handles.loadFilename_edit,'Visible','off');
if exist([filename,'.mat'])
    load(filename);
    handles.save.seiswav=seiswav;
    handles.save.seisatt=seisatt;
    handles.save.channelinfo=channelinfo;
    handles.AllFiles{1}=filename;
    set(handles.selectFiles,'String',filename);
    % Update handles structure
    guidata(hObject, handles);

    set(handles.list_button,'Enable','on');
    set(handles.map_button,'Enable','on');
    set(handles.plot_button,'Enable','on');

    set(handles.seis_no_end,'String',num2str(length(seiswav)));
    isparamerr(handles.seis_no_start,1,num2str(length(seiswav)),handles);
    set(handles.message,'String','File-Loading done');
else
    set(handles.message,'String','The given file not exists');
end
