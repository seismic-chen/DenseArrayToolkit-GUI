function varargout = SeisLister(varargin)
% SeisLister: a GUI program to list retrieved seismograms on IRIS and let user
% select interesting seismograms.
%
%   selection = SeisLister(mySeis, event) needs two inputs and returns one output. The
%   first input is mySeis, a retrieved MatSeismogram object. It should contain retrieved 
%   seismograms.The second input is a selected event, a structural array with the 
%   following stucture:
%           time: a string for event originTime, e.g. '2001-01-09T16:49:28.000Z'
%           longitude: a numerical number for event longitude, e.g. 167.1700
%           latitude: a numerical number for event latitude, e.g. -14.9280
%           depth: a numerical number for event depth, e.g. 103

%   The output, selection saves index numbers of seismograms user selected.
%
%   Also check MatSeismogram

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Aug, 2003


% Edit the above text to modify the response to help SeisLister

% Last Modified by GUIDE v2.5 25-Oct-2003 18:30:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeisLister_OpeningFcn, ...
                   'gui_OutputFcn',  @SeisLister_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before SeisLister is made visible.
function SeisLister_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for SeisLister
handles.output = hObject;
if nargin <=3
    msgbox('Need one retrieved MatSeismogram object...', 'Warning');
    delete(handles.figure1);
end

handles.mySeis=varargin{1};
handles.event=varargin{2};
handles.selection=[];
% Update handles structure
guidata(hObject, handles); 

if ~isempty(handles.mySeis)
    load_listbox(handles);
else
    beep;
    msgbox('Retrieve seismograms first...','Warning');
    return;
end
% UIWAIT makes SeisLister wait for user response (see UIRESUME)
guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeisLister_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles.selection)
    handles.selection = [];
else
    for i=1:length(handles.selection)
        handles.selection(i)=handles.selection(i)-1;
        if handles.selection(1) == 0 && i>1
            selection(i-1)=handles.selection(i);
        end
    end
    if handles.selection(1) == 0
        handles.selection = selection;
        clear selection;
    end
    
end
varargout{1} = handles.selection;
delete(handles.figure1);

% --- Executes on button press in station.
function station_Callback(hObject, eventdata, handles)
% hObject    handle to station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of station
load_listbox(handles);

% --- Executes on button press in lon.
function lon_Callback(hObject, eventdata, handles)
% hObject    handle to lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lon
load_listbox(handles);

% --- Executes on button press in lat.
function lat_Callback(hObject, eventdata, handles)
% hObject    handle to lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lat
load_listbox(handles);

% --- Executes on button press in dist.
function dist_Callback(hObject, eventdata, handles)
% hObject    handle to dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dist
load_listbox(handles);

% --- Executes on button press in azim.
function azim_Callback(hObject, eventdata, handles)
% hObject    handle to azim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of azim
load_listbox(handles);

% --- Executes on button press in starttime.
function starttime_Callback(hObject, eventdata, handles)
% hObject    handle to starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of starttime
load_listbox(handles);

% --- Executes on button press in endtime.
function endtime_Callback(hObject, eventdata, handles)
% hObject    handle to endtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of endtime
load_listbox(handles);

% --- Executes on button press in database.
function database_Callback(hObject, eventdata, handles)
% hObject    handle to database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of database
load_listbox(handles);

% --- Executes on button press in return_button.
function return_button_Callback(hObject, eventdata, handles)
% hObject    handle to return_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in site.
function site_Callback(hObject, eventdata, handles)
% hObject    handle to site (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of site
load_listbox(handles);

% --- Executes during object creation, after setting all properties.
function event_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in event_listbox.
function event_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to event_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns event_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from event_listbox
handles.selection = get(hObject, 'Value');
guidata(hObject, handles);

function load_listbox(handles);
% handles    structure with handles and user data (see GUIDATA)
% handles.mySeis    retrieved one MatSeismogram object, here it's mySeis

% check if seismograms are retrieved or not
mySeis=handles.mySeis;

% initialize fields
fields='';

attr=mySeis.seisatt;
% get attributions of retrieved seismograms
chan=mySeis.channelinfo;
% get channel information from retrieved MatSeismogram object
event=handles.event;

set(handles.display_Event, 'String', ['Event: originTime ', event.time,...
        '   Longitude ', num2str(event.longitude), '   Latitude ',...
        num2str(event.latitude),'   Depth ', num2str(event.depth)]);
% display event information

for i=1:length(attr)
    stringlistbox = '';     % initialize string
    if get(handles.network, 'Value') == get(handles.network, 'Max')
        stringlistbox = addspace(stringlistbox, attr(i).network,9);
        if i==1
            fields=addspace(fields,'Network',9);
        end
    end

    if get(handles.station, 'Value') == get(handles.station, 'Max')
        stringlistbox = addspace(stringlistbox, attr(i).station,9);
        if i==1
            fields=addspace(fields, 'Station',9);
        end
    end

    if get(handles.site, 'Value') == get(handles.site, 'Max')
        stringlistbox = addspace(stringlistbox, attr(i).site,6);
        if i==1
            fields=addspace(fields, 'Site',6);
        end
    end

    if get(handles.channel, 'Value') == get(handles.channel, 'Max')
        stringlistbox = addspace(stringlistbox, attr(i).channel,9);
        if i==1
            fields=addspace(fields, 'Channel',9);
        end
    end

    if get(handles.lon, 'Value') == get(handles.lon, 'Max')
       stringlistbox = addspace(stringlistbox, sprintf('%8.3f',chan(i).location.longitude),10);
        if i==1
            fields=addspace(fields, 'Longitude',10);
        end
    end
        
    if get(handles.lat, 'Value') == get(handles.lat, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%7.3f',chan(i).location.latitude), 10);
        if i==1
            fields=addspace(fields, 'Latitude',10);
        end
    end

    [delta,azim,backazim]=delaz(event.latitude,event.longitude, ...
            chan(i).location.latitude,  chan(i).location.longitude, 0);
    
    if get(handles.dist, 'Value') == get(handles.dist, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%3d',round(delta)),6);
        if i==1
            fields=addspace(fields, 'Dist',6);
        end
    end
        
    if get(handles.azim, 'Value') == get(handles.azim, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%3d',round(azim)),6);
        if i==1
            fields=addspace(fields, 'Azim',6);
        end
    end

    stx=timedif(event.time,attr(i).beginTime);
    etx= stx + attr(i).nPoint*attr(i).sampleIntv*0.001;
    if get(handles.starttime, 'Value') == get(handles.starttime, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%9.2f',stx), 12);
        if i==1
            fields=addspace(fields, '  startTime',12);
        end
    end
        
    if get(handles.endtime, 'Value') == get(handles.endtime, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%9.2f',etx),12);
        if i==1
            fields=addspace(fields, '  endTime',12);
        end
    end
       
    if get(handles.database, 'Value') == get(handles.database, 'Max')
        stringlistbox = addspace(stringlistbox, attr(i).from,8);
        if i==1
            fields=addspace(fields, 'Database',8);
        end
    end
    stringslistbox(i+1)=cellstr(stringlistbox);
end
stringslistbox(1)=cellstr(fields);
stringslistbox=stringslistbox';
set(handles.event_listbox, 'FontName', 'FixedWidth')
set(handles.event_listbox, 'FontSize', 12)
set(handles.event_listbox, 'String', stringslistbox);

function stringout = addspace(s1,s2,m)
% s1        the string to be appended
% s2        the string need appending
% m         the length of string s2 after padding with space
%           m >= length(s2)
% stringout = s1 + s2
lens2 = length(s2);
if lens2 < m
    for i=1:(m-lens2)
        s2 = [s2, ' '];
    end
end
stringout = [s1, s2];


% --- Executes on button press in network.
function network_Callback(hObject, eventdata, handles)
% hObject    handle to network (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of network
load_listbox(handles);

% --- Executes on button press in channel.
function channel_Callback(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of channel
load_listbox(handles);

