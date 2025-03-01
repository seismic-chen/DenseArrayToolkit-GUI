function varargout = SiteattRetriever(varargin)
% SITEATTRETRIEVER M-file for SiteattRetriever.fig
%      SITEATTRETRIEVER, by itself, creates a new SITEATTRETRIEVER or raises the existing
%      singleton*.
%
%      H = SITEATTRETRIEVER returns the handle to a new SITEATTRETRIEVER or the handle to
%      the existing singleton*.
%
%      SITEATTRETRIEVER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SITEATTRETRIEVER.M with the given input arguments.
%
%      SITEATTRETRIEVER('Property','Value',...) creates a new SITEATTRETRIEVER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SiteattRetriever_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SiteattRetriever_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SiteattRetriever

% Last Modified by GUIDE v2.5 01-Mar-2004 15:54:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SiteattRetriever_OpeningFcn, ...
                   'gui_OutputFcn',  @SiteattRetriever_OutputFcn, ...
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


% --- Executes just before SiteattRetriever is made visible.
function SiteattRetriever_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SiteattRetriever (see VARARGIN)
handles.output = hObject;
if nargin == 4
    handles.seisatt=varargin{1};
    handles.output=1;
else
    msgbox('Seisatt is required...','Warning');
    delete(handles.figure1);
    handles.output=0;
    return;
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SiteattRetriever wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SiteattRetriever_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);


% --- Executes during object creation, after setting all properties.
function min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_edit as text
%        str2double(get(hObject,'String')) returns contents of min_edit as a double
tmp=str2double(get(hObject, 'String'));
if isnan(tmp)
    beep; set(hObject, 'String', 'Error: Must be numeric!');
    set(hObject, 'ForegroundColor','red');
else
    set(hObject, 'ForegroundColor','black');
end

% --- Executes during object creation, after setting all properties.
function max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_edit as text
%        str2double(get(hObject,'String')) returns contents of max_edit as a double
tmp=str2double(get(hObject, 'String'));
if isnan(tmp)
    beep; set(hObject, 'String', 'Error: Must be numeric!');
    set(hObject, 'ForegroundColor','red');
else
    set(hObject, 'ForegroundColor','black');
end

% --- Executes during object creation, after setting all properties.
function num_freq_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function num_freq_edit_Callback(hObject, eventdata, handles)
% hObject    handle to num_freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_freq_edit as text
%        str2double(get(hObject,'String')) returns contents of num_freq_edit as a double
tmp=str2double(get(hObject, 'String'));
if isnan(tmp)
    beep; set(hObject, 'String', 'Error: Must be numeric!');
    set(hObject, 'ForegroundColor','red');
else
    set(hObject, 'ForegroundColor','black');
end

% --- Executes during object creation, after setting all properties.
function dir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dir_edit as text
%        str2double(get(hObject,'String')) returns contents of dir_edit as a double


% --- Executes during object creation, after setting all properties.
function property_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to property_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function property_edit_Callback(hObject, eventdata, handles)
% hObject    handle to property_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of property_edit as text
%        str2double(get(hObject,'String')) returns contents of property_edit as a double


% --- Executes on button press in disp_radio.
function disp_radio_Callback(hObject, eventdata, handles)
% hObject    handle to disp_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of disp_radio
if get(handles.disp_radio,'Value')==get(handles.disp_radio,'Max')
    set(handles.acc_radio,'Value',0.0);
    set(handles.default_radio,'Value',0.0);
    set(handles.vel_radio,'Value',0.0);
end

% --- Executes on button press in vel_radio.
function vel_radio_Callback(hObject, eventdata, handles)
% hObject    handle to vel_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vel_radio
if get(handles.vel_radio,'Value')==get(handles.vel_radio,'Max')
    set(handles.acc_radio,'Value',0.0);
    set(handles.default_radio,'Value',0.0);
    set(handles.disp_radio,'Value',0.0);
end

% --- Executes on button press in acc_radio.
function acc_radio_Callback(hObject, eventdata, handles)
% hObject    handle to acc_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of acc_radio
if get(handles.acc_radio,'Value')==get(handles.acc_radio,'Max')
    set(handles.disp_radio,'Value',get(handles.disp_radio,'Min'));
    set(handles.default_radio,'Value',get(handles.default_radio,'Min'));
    set(handles.vel_radio,'Value',get(handles.vel_radio,'Min'));
end

% --- Executes on button press in default_radio.
function default_radio_Callback(hObject, eventdata, handles)
% hObject    handle to default_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of default_radio
if get(handles.default_radio,'Value')==get(handles.default_radio,'Max')
    set(handles.acc_radio,'Value',0.0);
    set(handles.disp_radio,'Value',0.0);
    set(handles.vel_radio,'Value',0.0);
end

% --- Executes on button press in retrieve_button.
function retrieve_button_Callback(hObject, eventdata, handles)
% hObject    handle to retrieve_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
min=get(handles.min_edit,'String');
max=get(handles.max_edit,'String');
nfreq=get(handles.num_freq_edit,'String');
outdir=get(handles.dir_edit,'String');
prop=get(handles.property_edit,'String');
prop=which(prop);
if get(handles.disp_radio,'Value')==get(handles.disp_radio,'Max')
    unit='dis';
elseif get(handles.vel_radio,'Value')==get(handles.vel_radio,'Max')
    unit='vel';
elseif get(handles.acc_radio,'Value')==get(handles.acc_radio,'Max')
    unit='acc';
elseif get(handles.default_radio,'Value')==get(handles.default_radio,'Max')
    unit='def';
end

seisatt=handles.seisatt;
station=seisatt(1).station;
chan=seisatt(1).channel;
loc=seisatt(1).site;
net=seisatt(1).network;
nTotal=length(seisatt);
ab=com.isti.jevalresp.Run;
arg=[];

startyear='10000';
startdoy='10000';
endyear='0';
enddoy='0';
for i=1:nTotal
    beginTime=seisatt(i).beginTime;
    [d,t]=strtok(beginTime,'T');
    t=t(2:9);
    y=d(1:4);
    m=d(6:7);
    dy=num2str(doy(d));
    if str2num(y)<=str2num(startyear)
        startyear=y;
        if str2num(dy)<str2num(startdoy)
            startdoy=dy;
        end
    end
    if str2num(y)>=str2num(endyear)
        endyear=y;
        if str2num(dy)>str2num(enddoy)
            enddoy=dy;
        end
    end
    if isempty(strfind(net,seisatt(i).network))
        net=[net,',',seisatt(i).network];
    end
    % get the list of networks
    if isempty(strfind(station,seisatt(i).station))
        station=[station,',',seisatt(i).station];
    end
    % get the list of stations
    if isempty(strfind(chan,seisatt(i).channel))
        chan=[chan,',',seisatt(i).channel];
    end
    % get the list of channels
    if isempty(strfind(loc,seisatt(i).site))
        loc=[loc,',',seisatt(i).site];
    end
    % get the list of sites
end
 
 arg{1}=station;
 arg{2}=chan;
 arg{3}=startyear;
 arg{4}=startdoy;
 arg{5}=min;
 arg{6}=max;
 arg{7}=nfreq;
 arg{8}='-u';
 arg{9}=unit;
 arg{10}='-o';
 arg{11}=outdir;
 arg{12}='-v';
 arg{13}='-p';
 arg{14}=prop;
 arg{15}='-t';
 arg{16}='00:00:00';
 arg{17}='-n';
 arg{18}=net;
 arg{19}='-l';
 arg{20}=loc;
 arg{21}='-ey';
 arg{22}=endyear;
 arg{23}='-ed';
 arg{24}=enddoy;
 arg{25}='-et';
 arg{26}='23:59:59';
 clear seisatt;
if ~exist(prop)
    msgbox('DMC property file not exists in working Dir','Warning')
else
    ab.processAndOutput(arg);
end
uiresume(handles.figure1)
