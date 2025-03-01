function varargout = EventLister(varargin)
% EVENTLister: a GUI program to list found events on IRIS and let user
% select interesting events.
%
%   selection = event_Find(myEvent) needs one input and 
%   returns one output. The only input is myEvent, a MatEvent object.
%   It should contain retrieved events. The output, selection saves
%   index numbers of events user selected.
%
%   Also check MatEvent

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Aug, 2003


% Edit the above text to modify the response to help EventLister

% Last Modified by GUIDE v2.5 17-Jul-2003 15:30:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EventLister_OpeningFcn, ...
                   'gui_OutputFcn',  @EventLister_OutputFcn, ...
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


% --- Executes just before EventLister is made visible.
function EventLister_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for EventLister
handles.output = hObject;
if nargin <=3
    msgbox('Need one retrieved MatEvent object...', 'Warning');
    delete(handles.figure1);
end

handles.myEvent=varargin{1};
handles.selection=[];
% Update handles structure
guidata(hObject, handles); 

if ~isa(handles.myEvent, 'MatEvent')
    beep;
    msgbox('Input must be one object of MatEvent class','Warning');
    return;
end;

% load events which will be displayed in listbox
if ~isempty(handles.myEvent)
    load_listbox(handles);
else
    beep;
    msgbox('Retrieve event first...','Warning');
    return;
end
% UIWAIT makes EventLister wait for user response (see UIRESUME)
guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EventLister_OutputFcn(hObject, eventdata, handles)
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

% --- Executes on button press in name_checkbox.
function name_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to name_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of name_checkbox
load_listbox(handles);

% --- Executes on button press in longitude_checkbox.
function longitude_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to longitude_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of longitude_checkbox
load_listbox(handles);

% --- Executes on button press in latitude_checkbox.
function latitude_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to latitude_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of latitude_checkbox
load_listbox(handles);

% --- Executes on button press in depth_checkbox.
function depth_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to depth_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of depth_checkbox
load_listbox(handles);

% --- Executes on button press in magnitude_checkbox.
function magnitude_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to magnitude_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of magnitude_checkbox
load_listbox(handles);

% --- Executes on button press in origin_checkbox.
function origin_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to origin_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of origin_checkbox
load_listbox(handles);

% --- Executes on button press in catalog_checkbox.
function catalog_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to catalog_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of catalog_checkbox
load_listbox(handles);

% --- Executes on button press in contributor_checkbox.
function contributor_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to contributor_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of contributor_checkbox
load_listbox(handles);

% --- Executes on button press in update_pushbutton.
function update_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to update_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_listbox(handles);

% --- Executes on button press in return_button.
function return_button_Callback(hObject, eventdata, handles)
% hObject    handle to return_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in type_checkbox.
function type_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to type_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of type_checkbox
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
% handles.myEvent    retrieved one MatEvent object, here it's myEvent

% check if events are retrieved or not
if handles.myEvent.nEvent==0
    msgbox('Event information not retrieved yet, run Retrieve first...', 'Warning');
    delete(handles.figure1);
end

% initialize fields
fields='';

for i=1:handles.myEvent.nEvent
    event=handles.myEvent.preferredOrigin(i);   % get event information
    stringlistbox = '';     % initialize string
    if get(handles.no_checkbox, 'Value') == get(handles.no_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, num2str(i),4);
        if i==1
            fields=addspace(fields,'No',4);
        end
    end

    if get(handles.name_checkbox, 'Value') == get(handles.name_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, event.id,10);
        if i==1
            fields=addspace(fields, 'Name',10);
        end
    end
        
    if get(handles.longitude_checkbox, 'Value') == get(handles.longitude_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%8.3f',event.location.longitude),10);
        if i==1
            fields=addspace(fields, 'Longitude',10);
        end
    end
        
    if get(handles.latitude_checkbox, 'Value') == get(handles.latitude_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%7.3f',event.location.latitude), 10);
        if i==1
            fields=addspace(fields, 'Latitude',10);
        end
    end
        
    if get(handles.depth_checkbox, 'Value') == get(handles.depth_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, sprintf('%5.1f',event.location.depth),7);
        if i==1
            fields=addspace(fields, 'Depth',7);
        end
    end
        
    if get(handles.origin_checkbox, 'Value') == get(handles.origin_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, event.originTime, 27);
        if i==1
            fields=addspace(fields, 'originTime',27);
        end
    end
        
    if get(handles.catalog_checkbox, 'Value') == get(handles.catalog_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, event.catalog,8);
        if i==1
            fields=addspace(fields, 'Catalog',8);
        end
    end
       
    if get(handles.contributor_checkbox, 'Value') == get(handles.contributor_checkbox, 'Max')
        stringlistbox = addspace(stringlistbox, event.contributor,8);
        if i==1
            fields=addspace(fields, 'Contrb',8);
        end
    end

    if get(handles.magnitude_checkbox, 'Value') == get(handles.magnitude_checkbox, 'Max')
        for j=1:length(event.magnitude.value)
            stringlistbox = addspace(stringlistbox, sprintf('%4.1f',event.magnitude.value(j)),5);
        end
        if i==1
            fields=addspace(fields, ' Mag',5*length(event.magnitude.value));
        end
    end

    if get(handles.type_checkbox, 'Value') == get(handles.type_checkbox, 'Max')
        for j=1:length(event.magnitude.value)
            stringlistbox = addspace(stringlistbox, char(event.magnitude.type(j)),4);
        end
        if i==1
            fields=addspace(fields, 'Type',4*length(event.magnitude.value));
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


% --- Executes on button press in no_checkbox.
function no_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to no_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_checkbox
load_listbox(handles);

