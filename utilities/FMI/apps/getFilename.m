function varargout = getFilename(varargin)
% GETFILENAME M-file for getFilename.fig
%      GETFILENAME, by itself, creates a new GETFILENAME or raises the existing
%      singleton*.
%
%      H = GETFILENAME returns the handle to a new GETFILENAME or the handle to
%      the existing singleton*.
%
%      GETFILENAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETFILENAME.M with the given input arguments.
%
%      GETFILENAME('Property','Value',...) creates a new GETFILENAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getFilename_OpeningFunction gets called.  An
%      unrecognized property nameedit or invalid value makes property application
%      stop.  All inputs are passed to getFilename_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getFilename

% Last Modified by GUIDE v2.5 02-Feb-2004 15:40:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @getFilename_OpeningFcn, ...
                   'gui_OutputFcn',  @getFilename_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin>1 & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before getFilename is made visible.
function getFilename_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getFilename (see VARARGIN)

% Choose default command line output for getFilename
handles.output = hObject;
if length(varargin{1})>19
    handles.name=varargin{1}(1:19);
else
    handles.name=varargin{1};
end
handles.number=0;
set(handles.nameEdit,'String',handles.name);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes getFilename wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = getFilename_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.name;
varargout{2} = handles.number;
delete(handles.figure1);

% --- Executes during object creation, after setting all properties.
function nameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function nameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nameEdit as text
%        str2double(get(hObject,'String')) returns contents of nameEdit as a double


% --- Executes on button press in fullRetrieve.
function fullRetrieve_Callback(hObject, eventdata, handles)
% hObject    handle to fullRetrieve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fullRetrieve
if get(hObject,'Value')==1.0
    set(handles.numberEdit,'String','0');
    handles.number=0;
else
    set(handles.numberEdit,'String','100');
    handles.number=100;
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function numberEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numberEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function numberEdit_Callback(hObject, eventdata, handles)
% hObject    handle to numberEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numberEdit as text
%        str2double(get(hObject,'String')) returns contents of numberEdit as a double
set(handles.fullRetrieve,'Value',0.0);
number=str2num(get(hObject,'String'));
if number<1
    set(hObject,'String','1');
else
    number=round(number);
    set(hObject,'String',num2str(number));
end

% --- Executes on button press in okExit.
function okExit_Callback(hObject, eventdata, handles)
% hObject    handle to okExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.name=get(handles.nameEdit,'String');
if get(handles.fullRetrieve,'Value')==0.0
    handles.number=str2num(get(handles.numberEdit,'String'));
else
    handles.number=0;
end
% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1); 

