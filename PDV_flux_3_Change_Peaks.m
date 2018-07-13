function varargout = PDV_flux_3_Change_Peaks(varargin)
% PDV_FLUX_3_CHANGE_PEAKS MATLAB code for PDV_flux_3_Change_Peaks.fig
%      PDV_FLUX_3_CHANGE_PEAKS, by itself, creates a new PDV_FLUX_3_CHANGE_PEAKS or raises the existing
%      singleton*.
%
%      H = PDV_FLUX_3_CHANGE_PEAKS returns the handle to a new PDV_FLUX_3_CHANGE_PEAKS or the handle to
%      the existing singleton*.
%
%      PDV_FLUX_3_CHANGE_PEAKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PDV_FLUX_3_CHANGE_PEAKS.M with the given input arguments.
%
%      PDV_FLUX_3_CHANGE_PEAKS('Property','Value',...) creates a new PDV_FLUX_3_CHANGE_PEAKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PDV_flux_3_Change_Peaks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PDV_flux_3_Change_Peaks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PDV_flux_3_Change_Peaks

% Last Modified by GUIDE v2.5 08-Jul-2018 17:17:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PDV_flux_3_Change_Peaks_OpeningFcn, ...
                   'gui_OutputFcn',  @PDV_flux_3_Change_Peaks_OutputFcn, ...
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


% --- Executes just before PDV_flux_3_Change_Peaks is made visible.
function PDV_flux_3_Change_Peaks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PDV_flux_3_Change_Peaks (see VARARGIN)

% Choose default command line output for PDV_flux_3_Change_Peaks
handles.output = hObject;
h = findobj('Name','PDV_flux_3');
if ~isempty(h)
    main_GUI_data = guidata(h);
end
handles.main_GUI_data = main_GUI_data;
handles.state = '';
handles.channel = 1;
set(handles.files_list,'String',[handles.main_GUI_data.fnames']);
set(handles.files_list,'Value',handles.main_GUI_data.list_idx);
plot_int_cp(hObject,eventdata,handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PDV_flux_3_Change_Peaks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PDV_flux_3_Change_Peaks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in add_pnt.
function add_pnt_Callback(hObject, eventdata, handles)
% hObject    handle to add_pnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles.main_GUI_data;
n = h.list_idx;i = handles.channel;
if i==4
    j=3;
else
    j=i;
end
time = h.time{n};
amp = h.amplitude{n}(:,i);
if h.peaks_bool==1
    set(hObject,'Enable','off');
    set(handles.dlt_pnt,'Enable','off');
    [x_add,~] = ginput(1);
    set(hObject,'Enable','on');
    set(handles.dlt_pnt,'Enable','on');
    [~,idx] = min(abs(x_add-time));x_add = time(idx);
    h.xyPeaks{n}{j} = sortrows([h.xyPeaks{n}{j};x_add,amp(idx)]);
    handles.main_GUI_data.xyPeaks{n}{j} = h.xyPeaks{n}{j};
    [new_time,new_velocity,new_displacement] = vel_recalc(hObject,eventdata,handles);
    handles.main_GUI_data.lineout_time{n} = new_time;
    handles.main_GUI_data.velocity{n} = new_velocity;
    plot_int_cp(hObject,eventdata,handles);
    [handles.main_GUI_data.pressure{n},handles.main_GUI_data.flux{n},handles.main_GUI_data.fluence{n}] = calc_derived_data(hObject,eventdata,handles);
    handles.main_GUI_data.displacement{n} = new_displacement;
    guidata(hObject,handles);
    plot_int_cp(hObject,eventdata,handles);
end


% --- Executes on button press in dlt_pnt.
function dlt_pnt_Callback(hObject, eventdata, handles)
% hObject    handle to dlt_pnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles.main_GUI_data;
n = h.list_idx;i = handles.channel;
if i==4
    j=3;
else
    j=i;
end
time = h.time{n};
amp = h.amplitude{n}(:,i);
if h.peaks_bool==1
    set(hObject,'Enable','off');
    set(handles.dlt_pnt,'Enable','off');
    [x_rm,~] = ginput(1);
    set(hObject,'Enable','on');
    set(handles.dlt_pnt,'Enable','on');
    xvals = h.xyPeaks{n}{j}(:,1);
    [~,idx] = min(abs(xvals-x_rm));
    h.xyPeaks{n}{j}(idx,:) = [];
    handles.main_GUI_data.xyPeaks{n}{j} = h.xyPeaks{n}{j};
    [new_time,new_velocity,new_displacement] = vel_recalc(hObject,eventdata,handles);
    handles.main_GUI_data.lineout_time{n} = new_time;
    handles.main_GUI_data.velocity{n} = new_velocity;
    plot_int_cp(hObject,eventdata,handles);
    [handles.main_GUI_data.pressure{n},handles.main_GUI_data.flux{n},handles.main_GUI_data.fluence{n}] = calc_derived_data(hObject,eventdata,handles);
    handles.main_GUI_data.displacement{n} = new_displacement;
    guidata(hObject,handles);
    
end

% --- Executes on selection change in channel_list.
function channel_list_Callback(hObject, eventdata, handles)
% hObject    handle to channel_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channel_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel_list
channel = get(hObject,'Value');
switch channel
    case 1
        handles.channel =1;
    case 2
        handles.channel = 2;
    case 3
        handles.channel = 4;
end
guidata(hObject,handles);
plot_int_cp(hObject,eventdata,handles);
% --- Executes during object creation, after setting all properties.
function channel_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in files_list.
function files_list_Callback(hObject, eventdata, handles)
% hObject    handle to files_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns files_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from files_list
handles.main_GUI_data.list_idx = get(hObject,'Value');
guidata(hObject,handles);
plot_int_cp(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function files_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Name','PDV_flux_3');
guidata(h,handles.main_GUI_data);
close(PDV_flux_3_Change_Peaks);

% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(PDV_flux_3_Change_Peaks);

%% plots interferogram
    function int = plot_int_cp(hObject,eventdata,handles)
        h=handles.main_GUI_data;  n=h.list_idx;
        name = strsplit(h.fnames{n},'_Ch');
        name = name{1};
        amp = h.amplitude{n};
        axes(handles.main_axes); hold off;
        i = handles.channel;
        if i==4
            j=3;
        else
            j=i;
        end
        plot(h.time{n},amp(:,i)); hold on;
        if (h.method ==2 && h.peaks_bool ==1)
            xvals = h.xyPeaks{n}{j}(:,1);yvals = h.xyPeaks{n}{j}(:,2);
            plot(xvals,yvals,'.k','markersize',24)
        end
        
        xlabel('time (ns)');ylabel('volts (V)'); title(sprintf('%s: %s',name,handles.state));
        xlim([h.time{n}(h.t0{n})-10,h.time{n}(h.t0{n})+100]);
        chan_txt = sprintf('Channel %d',j);
        if h.peaks_bool ==0
            legend(chan_txt);
        else
            legend(chan_txt,'Registered Peaks');
        end
%% Begin subroutine section
    function [new_time,new_velocity,new_displacement]=vel_recalc(hObject,eventdata,handles)
        h = handles.main_GUI_data;
        n = h.list_idx;
        x0 = h.t0{n}-4;
        xyPeaks = h.xyPeaks{n};
        time = h.time{n};
        k = 1;
        for i=1:4
            if i==3
                continue
            end
        idx_0 = length(xyPeaks{k}(xyPeaks{k}(:,1)<=time(x0),1));
        velocity{k}=0.3875./diff(xyPeaks{k}(idx_0+1:end,1));
        xPeaks{k}=xyPeaks{k}(idx_0+1:end,1);xPeaks{k}(length(velocity{k}))=[];
        k = k+1;
        end
        velocity0 = 0;
        XYMAT=sortrows([xPeaks{1},velocity{1};xPeaks{2},velocity{2};xPeaks{3},velocity{3};time(x0),velocity0]);
        new_time=XYMAT(:,1);
        new_velocity=XYMAT(:,2)/h.window_val;
        x = smooth(new_velocity(3:end),3);
        new_velocity(3:end) = x;
        
        %recalculate displacement
        displacement1=cell(length(xPeaks{1}),1);
        displacement2=cell(length(xPeaks{2}),1);
        displacement3=cell(length(xPeaks{3}),1);
        [~, minimum_index] = min(h.phase{n});
        time = h.time{n};
        for i=1:length(xPeaks{1}(:,1))
            if i==1
                displacement1{i}=(1.55/(4*h.window_val));
            elseif xPeaks{1}(i) < time(minimum_index)
                displacement1{i}=displacement1{i-1}+(1.55/(4*h.window_val));
            else
                displacement1{i}=displacement1{i-1}-(1.55/(4*h.window_val));
            end
        end
        for i=1:length(xPeaks{2})
            if i==1
                displacement2{i}=(1.55/(4*h.window_val));
            elseif xPeaks{2}(i) < time(minimum_index)
                displacement2{i}=displacement2{i-1}+(1.55/(4*h.window_val));
            else
                displacement2{i}=displacement2{i-1}-(1.55/(4*h.window_val));
            end
        end
        for i=1:length(xPeaks{3})
            if i==1
                displacement3{i}=(1.55/(4*h.window_val));
            elseif xPeaks{3}(i) < time(minimum_index)
                displacement3{i}=displacement3{i-1}+(1.55/(4*h.window_val));
            else
                displacement3{i}=displacement3{i-1}-(1.55/(4*h.window_val));
            end
        end
        displacement_c=cell(3,2);
        displacement_c{1,1}=xPeaks{1};displacement_c{2,1}=xPeaks{2};displacement_c{3,1}=xPeaks{3};
        displacement_c{1,2}=cell2mat(displacement1);displacement_c{2,2}=cell2mat(displacement2);displacement_c{3,2}=cell2mat(displacement3);
        displacement2c=cell2mat(displacement_c);
        displacement2c=sortrows(displacement2c);
        displacementTime = displacement2c(:,1);
        displacement=displacement2c(:,2);
        new_displacement=smooth(displacementTime,displacement,0.1,'rloess');
        new_displacement = [0;new_displacement];
        %% Recalculate pressure,flux and fluence
function [pressure,flux,fluence] = calc_derived_data(hObject,eventdata,handles)
        i=1;
        j=1;
        h = handles.main_GUI_data;
        n = h.list_idx;
        velocity_lineout_fit = h.velocity{n};
        time_lineout = h.lineout_time{n};
        fluence(j,1) = 0;
        pressure(j,1)=0;
        switch h.sample_ID
            case 1
            rho = 2.230;
        while i<length(velocity_lineout_fit)
            j=j+1;
            if abs(velocity_lineout_fit(i,1)) <= 0.568
                m = 1.861;%% note: these are specific to GLASS, from the glass hugoniot --> that is why the constants change for different velocities 
                b = 3.879;
                flux(j,1)=find_flux(m,b,rho,velocity_lineout_fit(i,1));% the flux is = force * dA , where A is area 
                pressure((j),1)=find_pressure(m,b,rho,velocity_lineout_fit((i),1));
                fluence(j,1)=fluence((j-1),1)+trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
            elseif abs(velocity_lineout_fit(i,1)) > 1.83
                m = 1.269;
                b = 2.3925;
                flux(j,1)=find_flux(m,b,rho,velocity_lineout_fit(i,1)); 
                pressure((j),1)=find_pressure(m,b,rho,velocity_lineout_fit((i),1));
                fluence(j,1)=fluence((j-1),1)+trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
            else 
                m = -0.175;
                b = 5.0344; 
                flux(j,1)=find_flux(m,b,rho,velocity_lineout_fit(i,1));% the flux is = force * dA , where A is area 
                pressure((j),1)=find_pressure(m,b,rho,velocity_lineout_fit((i),1));
                fluence(j,1)=fluence((j-1),1)+trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
            end 
            i=i+1;
        end
            case 2 %Nitromethane
                rho = 1.125; m = 1.64; b = 1.65;
                while i<length(velocity_lineout_fit)
                    j=j+1;
                    flux(j,1) = find_flux(m,b,rho,velocity_lineout_fit(i,1));
                    pressure(j,1) = find_pressure(m,b,rho,velocity_lineout_fit(i,1));
                    fluence(j,1) = fluence((j-1),1)+trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
                    i = i+1;
                end
            case 3 %Sapphire
                rho = 4.0; m = 0.957; b = 8.74;
                while i<length(velocity_lineout_fit)
                    j=j+1;
                    flux(j,1) = find_flux(m,b,rho,velocity_lineout_fit(i,1));
                    pressure(j,1) = find_pressure(m,b,rho,velocity_lineout_fit(i,1));
                    fluence(j,1) = fluence((j-1),1)+trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
                    i = i+1;
                end
            case 4 %CaF2
                %Double check this
                rho = 3.18; m = 1.18; b = 5.15;
                while i<length(velocity_lineout_fit)
                    j=j+1;
                    flux(j,1) = find_flux(m,b,rho,velocity_lineout_fit(i,1));
                    pressure(j,1) = find_pressure(m,b,rho,velocity_lineout_fit(i,1));
                    fluence(j,1) = fluence((j-1),1)+trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
                    i = i+1;
                end
        end
        
    function new_displacement = recalc_disp(hObject,eventdata,handles)
        h = handles.main_GUI_data;
        n = h.list_idx;
        xPeaks = h.xyPeaks{n};
        displacement1=cell(length(xPeaks{1}(:,1)),1);
        displacement2=cell(length(xPeaks{2}(:,1)),1);
        displacement3=cell(length(xPeaks{3}(:,1)),1);
        [~, minimum_index] = min(h.phase{n});
        time = h.time{n};
        for i=1:length(xPeaks{1}(:,1))
            if i==1
                displacement1{i}=(1.55/(4*h.window_val));
            elseif xPeaks{1}(i,1) < time(minimum_index)
                displacement1{i}=displacement1{i-1}+(1.55/(4*h.window_val));
            else
                displacement1{i}=displacement1{i-1}-(1.55/(4*h.window_val));
            end
        end
        for i=1:length(xPeaks{2}(:,1))
            if i==1
                displacement2{i}=(1.55/(4*h.window_val));
            elseif xPeaks{2}(i,1) < time(minimum_index)
                displacement2{i}=displacement2{i-1}+(1.55/(4*h.window_val));
            else
                displacement2{i}=displacement2{i-1}-(1.55/(4*h.window_val));
            end
        end
        for i=1:length(xPeaks{3}(:,1))
            if i==1
                displacement3{i}=(1.55/(4*h.window_val));
            elseif xPeaks{3}(i,1) < time(minimum_index)
                displacement3{i}=displacement3{i-1}+(1.55/(4*h.window_val));
            else
                displacement3{i}=displacement3{i-1}-(1.55/(4*h.window_val));
            end
        end
        displacement_c=cell(3,2);
        displacement_c{1,1}=xPeaks{1};displacement_c{2,1}=xPeaks{2};displacement_c{3,1}=xPeaks{3};
        displacement_c{1,2}=cell2mat(displacement1);displacement_c{2,2}=cell2mat(displacement2);displacement_c{3,2}=cell2mat(displacement3);
        displacement2c=cell2mat(displacement_c);
        displacement2c=sortrows(displacement2c);
        displacementTime = displacement2c(:,1);
        displacement=displacement2c(:,2);
        new_displacement=smooth(displacementTime,displacement,0.1,'rloess');
        new_displacement = [0;new_displacement];

%% Find flux and fluence equations;        
    function flux = find_flux(m,b,rho,velocity)
        flux = 0.5*rho*(b+m*velocity)*velocity^2; % the flux is = force * dA , where A is area


%%function to calculate pressure given m,p,b,and vel
    function pressure = find_pressure(m,b,rho,velocity)
        pressure = rho*(b+m*velocity)*velocity;