function varargout = PDV_flux_3_Settings(varargin)
% PDV_FLUX_3_SETTINGS MATLAB code for PDV_flux_3_Settings.fig
%      PDV_FLUX_3_SETTINGS, by itself, creates a new PDV_FLUX_3_SETTINGS or raises the existing
%      singleton*.
%
%      H = PDV_FLUX_3_SETTINGS returns the handle to a new PDV_FLUX_3_SETTINGS or the handle to
%      the existing singleton*.
%
%      PDV_FLUX_3_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PDV_FLUX_3_SETTINGS.M with the given input arguments.
%
%      PDV_FLUX_3_SETTINGS('Property','Value',...) creates a new PDV_FLUX_3_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PDV_flux_3_Settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PDV_flux_3_Settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PDV_flux_3_Settings

% Last Modified by GUIDE v2.5 12-Jul-2018 13:59:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PDV_flux_3_Settings_OpeningFcn, ...
                   'gui_OutputFcn',  @PDV_flux_3_Settings_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PDV_flux_3_Settings is made visible.
function PDV_flux_3_Settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PDV_flux_3_Settings (see VARARGIN)

% Choose default command line output for PDV_flux_3_Settings
handles.output = hObject;
h = findobj('Name','PDV_flux_3');
if ~isempty(h)
    main_GUI_data = guidata(h);
end
handles.main_GUI_data = main_GUI_data;
if handles.main_GUI_data.method ==1
    set(handles.STFT_button,'Value',1);
    set(handles.fringe_button,'Value',0);
else
    set(handles.STFT_button,'Value',0);
    set(handles.fringe_button,'Value',1);
end
switch handles.main_GUI_data.time_resolution
    case 2
        set(handles.sel_tw,'Value',1);
    case 5
        set(handles.sel_tw,'Value',2);
    case 10
        set(handles.sel_tw,'Value',3);
    case 15
        set(handles.sel_tw,'Value',4);
end
set(handles.f_sample,'String',string(handles.main_GUI_data.sample_rate));
set(handles.v_cutoff,'String',string(handles.main_GUI_data.velocity_cutoff_low));
set(handles.t_offset,'String',string(handles.main_GUI_data.scope_offset));
set(handles.window_menu,'Value',handles.main_GUI_data.window_material);
set(handles.sample_material,'Value',handles.main_GUI_data.sample_ID);
set(handles.threshold_edit,'String',handles.main_GUI_data.threshold_val);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PDV_flux_3_Settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PDV_flux_3_Settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function t_offset_Callback(hObject, eventdata, handles)
% hObject    handle to t_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_offset as text
%        str2double(get(hObject,'String')) returns contents of t_offset as a double
xi = handles.main_GUI_data.scope_offset;
handles.main_GUI_data.scope_offset = str2double(get(hObject,'String'));
diff = handles.main_GUI_data.scope_offset-xi;
try
handles.main_GUI_data.time = cellfun(@(x) x+diff, handles.main_GUI_data.time, 'UniformOutput', false);
catch
end
try
handles.main_GUI_data.lineout_time = cellfun(@(x) x+diff, handles.main_GUI_data.lineout_time, 'UniformOutput', false);
catch
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function t_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sel_tw.
function sel_tw_Callback(hObject, eventdata, handles)
% hObject    handle to sel_tw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sel_tw contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_tw
x = get(hObject,'Value');
switch x
    case 1
        handles.main_GUI_data.time_resolution = 2;
    case 2
        handles.main_GUI_data.time_resolution = 5;
    case 3
        handles.main_GUI_data.time_resolution = 10;
    case 4
        handles.main_GUI_data.time_resolution = 15;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sel_tw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_tw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_sample_Callback(hObject, eventdata, handles)
% hObject    handle to f_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_sample as text
%        str2double(get(hObject,'String')) returns contents of f_sample as a double
handles.main_GUI_data.sample_rate = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function f_sample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to v_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v_cutoff as text
%        str2double(get(hObject,'String')) returns contents of v_cutoff as a double
handles.main_GUI_data.velocity_cutoff_low = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function v_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in STFT_button.
function STFT_button_Callback(hObject, eventdata, handles)
% hObject    handle to STFT_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of STFT_button
if get(hObject,'Value')==1
    set(handles.fringe_button,'Value',0);
    handles.main_GUI_data.method = 1;
end
guidata(hObject,handles);


 
 
% --- Executes on button press in fringe_button.
function fringe_button_Callback(hObject, eventdata, handles)
% hObject    handle to fringe_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fringe_button
if get(hObject,'Value')==1
    set(handles.STFT_button,'Value',0);
    handles.main_GUI_data.method = 2;
end
guidata(hObject,handles);

% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%main_GUI=guihandles('Tag','Emission_spectrometer_GUI');
m = handles.main_GUI_data.method;
if m==1
    set(handles.main_GUI_data.method_text,'String',sprintf('Currently Plotted: \n STFT'));
elseif m==2
    set(handles.main_GUI_data.method_text,'String',sprintf('Currently Plotted: \n Fringe Analysis'));
end
h = findobj('Name','PDV_flux_3');
guidata(h,handles.main_GUI_data);
close(PDV_flux_3_Settings);


% --- Executes on selection change in window_menu.
function window_menu_Callback(hObject, eventdata, handles)
% hObject    handle to window_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns window_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from window_menu
handles.main_GUI_data.window_material = get(hObject,'Value');
switch get(hObject,'Value')
    case 1
        handles.main_GUI_data.window_val = 1.065;
    case 2
        handles.main_GUI_data.window_val = 1.7462;
    case 3
        handles.main_GUI_data.window_val = 1.25;
    case 4
        handles.main_GUI_data.window_val = 1.3827;
    case 5
        handles.main_GUI_data.window_val = 1.5278;
end
%Window material choice for PDV window correction. Changing these will
%allow the correction factors to be input into the main function
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function window_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sample_material.
function sample_material_Callback(hObject, eventdata, handles)
% hObject    handle to sample_material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sample_material contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sample_material
handles.main_GUI_data.sample_ID = get(hObject,'Value');
%Selects shocked material to properly evaluate the flux values from
%velocity history.
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sample_material_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_edit as text
%        str2double(get(hObject,'String')) returns contents of threshold_edit as a double
handles.main_GUI_data.threshold_val = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function threshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
