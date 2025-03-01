function varargout = event_Find(varargin)
% EVENT_FIND: a GUI program to search all available events on IRIS
%   given a set of parameters, such as earthquake location, magnitude,
%   depth, etc. It requires event_Find.fig
%
%   [param, event, selection]=event_Find(myEvent) needs one input and 
%   returns three outputs. The only input is myEvent, a MatEvent object. 
%   The first output, param, saves parameters one input for searching 
%   earthquakes. The second one, event, is a MatEvent object 
%   which save all found events. The third one is selection saving index 
%   numbers of events one selected.
%
%   Also check MatEvent

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Sep, 2003

% Edit the above text to modify the response to help event_Find

% Last Modified by GUIDE v2.5 29-Jul-2003 12:58:55

% Begin initialization code - DO NOT MAG_END
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @event_Find_OpeningFcn, ...
                   'gui_OutputFcn',  @event_Find_OutputFcn, ...
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
% End initialization code - DO NOT MAG_END

% --- Executes just before event_Find is made visible.
function event_Find_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to event_Find (see VARARGIN)

% Choose default command line output for event_Find
handles.myEvent=varargin{1};
handles.selection=[];
handles.p=[];

% Update handles structure
guidata(hObject, handles);

if nargin <=3 || ~isa(varargin{1}, 'MatEvent')
    beep;
    msgbox('Input argument must be one MatEvent object...', 'Warning');
    handles.myEvent=[];
    return
end

% UIWAIT makes event_Find wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = event_Find_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.p;
varargout{2} = handles.myEvent;
varargout{3} = handles.selection;
delete(handles.figure1);

% --- Executes during object creation, after setting all properties.
function lat_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: mag_end controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function lat_start_Callback(hObject, eventdata, handles)
% hObject    handle to lat_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lat_start as text
%        str2double(get(hObject,'String')) returns contents of lat_start as a double
isparamerr(hObject, -90, 90, handles);

% --- Executes during object creation, after setting all properties.

function long_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to long_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: mag_end controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function long_start_Callback(hObject, eventdata, handles)
% hObject    handle to long_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of long_start as text
%        str2double(get(hObject,'String')) returns contents of long_start as a double
% --- Check if input is correct
isparamerr(hObject, -180, 180, handles);

% --- Executes during object creation, after setting all properties.
function lat_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: mag_end controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function lat_end_Callback(hObject, eventdata, handles)
% hObject    handle to lat_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lat_end as text
%        str2double(get(hObject,'String')) returns contents of lat_end as a double
% --- Get the ending latitude used in calculating area
isparamerr(hObject, -90, 90, handles);
guidata(hObject, handles);

function time_start_year_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_start_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: mag_end controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function time_start_year_Callback(hObject, eventdata, handles)
% hObject    handle to time_start_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_start_year as text
%        str2double(get(hObject,'String')) returns contents of time_start_year as a double
%handles.timestart=get(hObject, 'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function mag_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mag_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: mag_end controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function mag_start_Callback(hObject, eventdata, handles)
% hObject    handle to mag_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mag_start as text
%        str2double(get(hObject,'String')) returns contents of mag_start as a double
% --- Check if input is correct
isparamerr(hObject, 0, 12, handles);

% --- Executes during object creation, after setting all properties.
function mag_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mag_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: mag_end controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function mag_end_Callback(hObject, eventdata, handles)
% hObject    handle to mag_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of mag_end as text
%        str2double(get(hObject,'String')) returns contents of mag_end as a double
% --- Check if input is correct
isparamerr(hObject, 0, 12, handles);

% --- Executes during object creation, after setting all properties.
function event_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: mag_end controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function event_type_Callback(hObject, eventdata, handles)
% hObject    handle to event_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_type as text
%        str2double(get(hObject,'String')) returns contents of event_type as a double
% --- Get the event type


% --- Executes on button press in search_button.
function search_button_Callback(hObject, eventdata, handles)
% hObject    handle to search_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message_edit, 'String','Searching available events by given arguments...');

% --- Get all parameters
latstart = get(handles.lat_start,'String');
latend = get(handles.lat_end,'String');
longstart = get(handles.long_start,'String');
longend = get(handles.long_end,'String');
handles.p.area=[latstart, ' ', latend, ' ', longstart, ' ', longend];
handles.p.area=[str2num(handles.p.area)];
guidata(hObject, handles);

% get magnitude range
magstart = get(handles.mag_start,'String');
magend = get(handles.mag_end,'String');
handles.p.mag=[magstart, ' ', magend];
handles.p.mag=[str2num(handles.p.mag)];
guidata(hObject, handles);

% get depth range
depthstart = get(handles.depth_start,'String');
depthend = get(handles.depth_end,'String');
handles.p.depth=[depthstart, ' ', depthend];
handles.p.depth=[str2num(handles.p.depth)];
guidata(hObject, handles);

% get time start
timestartyear = addpad(get(handles.time_start_year,'String'),4);
timestartmon = addpad(get(handles.time_start_mon,'String'),2);
timestartday = addpad(get(handles.time_start_day,'String'),2);
timestarthour = addpad(get(handles.time_start_hour,'String'),2);
timestartmin = addpad(get(handles.time_start_min,'String'),2);
timestartsec = addpad(get(handles.time_start_sec,'String'),3);
timestart = [timestartyear, '-', timestartmon, '-', timestartday, 'T', timestarthour, ':', timestartmin, ':', timestartsec, 'Z'];

% get time end
timeendyear = addpad(get(handles.time_end_year,'String'),4);
timeendmon = addpad(get(handles.time_end_mon,'String'),2);
timeendday = addpad(get(handles.time_end_day,'String'),2);
timeendhour = addpad(get(handles.time_end_hour,'String'),2);
timeendmin = addpad(get(handles.time_end_min,'String'),2);
timeendsec = addpad(get(handles.time_end_sec,'String'),3);
timeend = [timeendyear, '-', timeendmon, '-', timeendday, 'T', timeendhour, ':', timeendmin, ':', timeendsec, 'Z'];
handles.p.time={timestart, timeend};
guidata(hObject, handles);

% get event type
i=1;
ptype=get(handles.event_type,'String');
while (any(ptype))
  [chopped,ptype] = strtok(ptype);
% here, delimiter is whitespace
  handles.p.type{i} = upper(chopped);
  i=i+1;
end
guidata(hObject, handles);

% get catalog
handles.p.catalog={upper(get(handles.event_catalog,'String'))};
guidata(hObject, handles);

% set buttons' properties
set(handles.search_button, 'Enable', 'off');
set(handles.list_button, 'Enable', 'off');
set(handles.map_button, 'Enable', 'off');
set(handles.retrieve_button, 'Enable', 'off');

handles.myEvent=find_event(handles.myEvent,'Area',handles.p.area,'Depth', handles.p.depth, 'magnitude',handles.p.mag, 'searchTypes',handles.p.type,'timeRange',handles.p.time, 'catalogs', handles.p.catalog);
set(handles.message_edit, 'String', [num2str(handles.myEvent.nEvent), ' events found']);
guidata(hObject, handles);
% activate retrieve button
set(handles.retrieve_button, 'Enable', 'on');
set(handles.search_button, 'Enable', 'on');

% --- Executes on button press in list_button.
function list_button_Callback(hObject, eventdata, handles)
% hObject    handle to list_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message_edit, 'String','Listing events...');
selection = EventLister(handles.myEvent);
if isempty(selection)
    selection=[1:handles.myEvent.nEvent];
end
handles.selection=selection;
guidata(hObject, handles);
set(handles.message_edit, 'String','Events-Listing done');

% --- Executes on button press in map_button.
function map_button_Callback(hObject, eventdata, handles)
% hObject    handle to map_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message_edit, 'String',' Mapping events...');
EventMaper(handles.myEvent, handles.p.area, handles.selection);
set(handles.message_edit, 'String','Event-Mapping done');

% --- Executes on button press in retrieve_button.
function retrieve_button_Callback(hObject, eventdata, handles)
% hObject    handle to retrieve_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.message_edit, 'String','Retrieving events...');
handles.myEvent=retrieve_event_info(handles.myEvent);
handles.selection=[1:handles.myEvent.nEvent];
guidata(hObject, handles);
set(handles.list_button, 'Enable','on');
set(handles.map_button,'Enable','on');
set(handles.message_edit, 'String','Events-retrieving done');

% --- Executes during object creation, after setting all properties.
function long_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to long_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function long_end_Callback(hObject, eventdata, handles)
% hObject    handle to long_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of long_end as text
%        str2double(get(hObject,'String')) returns contents of long_end as a double
% --- Check if input is right
isparamerr(hObject, -180, 180, handles);


% --- Executes during object creation, after setting all properties.
function message_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to message_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function message_edit_Callback(hObject, eventdata, handles)
% hObject    handle to message_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of message_edit as text
%        str2double(get(hObject,'String')) returns contents of message_edit as a double


% --- Executes during object creation, after setting all properties.
function time_start_mon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_start_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_start_mon_Callback(hObject, eventdata, handles)
% hObject    handle to time_start_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_start_mon as text
%        str2double(get(hObject,'String')) returns contents of time_start_mon as a double
% --- Check if input is correct
isparamerr(hObject, 0, 12, handles);


% --- Executes during object creation, after setting all properties.
function time_start_day_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_start_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_start_day_Callback(hObject, eventdata, handles)
% hObject    handle to time_start_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_start_day as text
%        str2double(get(hObject,'String')) returns contents of time_start_day as a double
% --- Check if input is correct
isparamerr(hObject, 0, 31, handles);


% --- Executes during object creation, after setting all properties.
function time_start_hour_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_start_hour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_start_hour_Callback(hObject, eventdata, handles)
% hObject    handle to time_start_hour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_start_hour as text
%        str2double(get(hObject,'String')) returns contents of time_start_hour as a double
% --- Check if input is correct
isparamerr(hObject, 0, 24, handles);


% --- Executes during object creation, after setting all properties.
function time_start_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_start_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_start_min_Callback(hObject, eventdata, handles)
% hObject    handle to time_start_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_start_min as text
%        str2double(get(hObject,'String')) returns contents of time_start_min as a double
% --- Check if input is correct
isparamerr(hObject, 0, 60, handles);


% --- Executes during object creation, after setting all properties.
function time_start_sec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_start_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_start_sec_Callback(hObject, eventdata, handles)
% hObject    handle to time_start_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_start_sec as text
%        str2double(get(hObject,'String')) returns contents of time_start_sec as a double
% --- Check if input is correct
isparamerr(hObject, 0, 60, handles);


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_start_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to time_start_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_start_year as text
%        str2double(get(hObject,'String')) returns contents of time_start_year as a double


% --- Executes during object creation, after setting all properties.
function time_end_mon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_end_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_end_mon_Callback(hObject, eventdata, handles)
% hObject    handle to time_end_mon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_end_mon as text
%        str2double(get(hObject,'String')) returns contents of time_end_mon as a double
% --- Check if input is correct
isparamerr(hObject, 0, 12, handles);


% --- Executes during object creation, after setting all properties.
function time_end_day_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_end_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_end_day_Callback(hObject, eventdata, handles)
% hObject    handle to time_end_day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_end_day as text
%        str2double(get(hObject,'String')) returns contents of time_end_day as a double
% --- Check if input is correct
isparamerr(hObject, 0, 31, handles);


% --- Executes during object creation, after setting all properties.
function time_end_hour_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_end_hour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_end_hour_Callback(hObject, eventdata, handles)
% hObject    handle to time_end_hour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_end_hour as text
%        str2double(get(hObject,'String')) returns contents of time_end_hour as a double
% --- Check if input is correct
isparamerr(hObject, 0, 24, handles);

% --- Executes during object creation, after setting all properties.
function time_end_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_end_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_end_min_Callback(hObject, eventdata, handles)
% hObject    handle to time_end_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_end_min as text
%        str2double(get(hObject,'String')) returns contents of time_end_min as a double
% --- Check if input is correct
isparamerr(hObject, 0, 60, handles);


% --- Executes during object creation, after setting all properties.
function time_end_sec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_end_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_end_sec_Callback(hObject, eventdata, handles)
% hObject    handle to time_end_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_end_sec as text
%        str2double(get(hObject,'String')) returns contents of time_end_sec as a double
% --- Check if input is correct
isparamerr(hObject, 0, 60, handles);


% --- Executes during object creation, after setting all properties.
function time_end_year_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_end_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_end_year_Callback(hObject, eventdata, handles)
% hObject    handle to time_end_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_end_year as text
%        str2double(get(hObject,'String')) returns contents of time_end_year as a double


%------------------------------------------------------------------------
% This function is used to tell if your input parameter is correct in GUI
% The correct parameter should be numeric and have an upper bound number
% eg, month must be numeric and less than 12. In order to run this function
% one 'edit' component with a 'message_edit' handle is necessary to display 
% warning message.
function isparamerr(hObject, lmit, umit, handles)
% Usage:    bObject,    the handle of one GUI edit component
%           umit,       the number's upper limit
%           lmit,       the number's lower limit
%           handles,    the handle struction
tmp=str2double(get(hObject, 'String'));
if isnan(tmp)
    beep; set(handles.message_edit, 'String', 'Error: Must be numeric !');
    set(hObject, 'ForegroundColor','red');
else
    if tmp > umit 
        beep;
        set(hObject, 'String', num2str(umit));
        set(hObject, 'ForegroundColor','red');
        set(handles.message_edit, 'String', ['Must be less than ', num2str(umit), ' !']);
    elseif tmp <lmit
        beep;
        set(hObject, 'String', num2str(lmit));
        set(hObject, 'ForegroundColor','red');
        set(handles.message_edit, 'String', ['Must be greater than ', num2str(lmit), ' !']);
    else
        set(hObject, 'ForegroundColor','black');
        set(handles.message_edit, 'String', '');
    end
end
clear tmp;

%----------------------------------------------------------------
% function readparam(handles)
function readparam(hObject, handles)
latstart = get(handles.lat_start,'String');
latend = get(handles.lat_end,'String');
longstart = get(handles.long_start,'String');
longend = get(handles.long_end,'String');
handles.p.area=[latstart, ' ', latend, ' ', longstart, ' ', longend];
handles.p.area=[str2num(handles.p.area)];
% get magnitude range
magstart = get(handles.mag_start,'String');
magend = get(handles.mag_end,'String');
handles.p.mag=[magstart, ' ', magend];
handles.p.mag=[str2num(handles.p.mag)];

% get time start
timestartyear = get(handles.time_start_year,'String');
timestartmon = get(handles.time_start_mon,'String');
timestartday = get(handles.time_start_day,'String');
timestarthour = get(handles.time_start_hour,'String');
timestartmin = get(handles.time_start_min,'String');
timestartsec = get(handles.time_start_sec,'String');
timestart = [timestartyear, '-', timestartmon, '-', timestartday, 'T', timestarthour, ':', timestartmin, ':', timestartsec, 'Z'];
% get time end
timeendyear = get(handles.time_end_year,'String');
timeendmon = get(handles.time_end_mon,'String');
timeendday = get(handles.time_end_day,'String');
timeendhour = get(handles.time_end_hour,'String');
timeendmin = get(handles.time_end_min,'String');
timeendsec = get(handles.time_end_sec,'String');
timeend = [timeendyear, '-', timeendmon, '-', timeendday, 'T', timeendhour, ':', timeendmin, ':', timeendsec, 'Z'];
handles.p.time={timestart, timeend};
% get event type
handles.p.eventtype={get(handles.event_type,'String')};

% update handles structure
guidata(hObject, handles);

% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function depth_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depth_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function depth_start_Callback(hObject, eventdata, handles)
% hObject    handle to depth_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depth_start as text
%        str2double(get(hObject,'String')) returns contents of depth_start as a double
isparamerr(hObject, 0, 1000, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function depth_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depth_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function depth_end_Callback(hObject, eventdata, handles)
% hObject    handle to depth_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depth_end as text
%        str2double(get(hObject,'String')) returns contents of depth_end as a double
isparamerr(hObject, 0, 1000, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function event_catalog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function event_catalog_Callback(hObject, eventdata, handles)
% hObject    handle to event_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_catalog as text
%        str2double(get(hObject,'String')) returns contents of event_catalog as a double


function outstr=addpad(instr, lenstr)
% instr     source string
% lenstr    length of output string

if length(instr)<lenstr
    for i=1:lenstr-length(instr)
        instr=['0',instr];
    end
end
outstr=instr;