function varargout = PDV_flux_3(varargin)
% PDV_FLUX_3 MATLAB code for PDV_flux_3.fig
%      PDV_FLUX_3, by itself, creates a new PDV_FLUX_3 or raises the existing
%      singleton*.
%
%      H = PDV_FLUX_3 returns the handle to a new PDV_FLUX_3 or the handle to
%      the existing singleton*.
%
%   n (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PDV_flux_3

% Last Modified by GUIDE v2.5 28-Aug-2018 15:49:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PDV_flux_3_OpeningFcn, ...
                   'gui_OutputFcn',  @PDV_flux_3_OutputFcn, ...
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


% --- Executes just before PDV_flux_3 is made visible.
function PDV_flux_3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PDV_flux_3 (see VARARGIN)

% Choose default command line output for PDV_flux_3
handles.xpeaks = {};
handles.calc_bool = 0;
handles.peaks_bool = 0;
handles.output = hObject;
handles.scope_offset = 23;
handles.time_resolution = 5;
handles.sample_rate = 0.08;
handles.velocity_cutoff_low = 0;
handles.method = 1;
handles.sample_ID = 1;
handles.window_material = 1;
handles.window_val = 1.0627;
handles.threshold_val = 3;
handles.list_idx = 0;
handles.plot_idx = 1;
handles.save_timing_bool = 0;
handles.fnames = {};
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PDV_flux_3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PDV_flux_3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% --------------------------------------------------------------------
function Load_data_tag_Callback(hObject, eventdata, handles)
% hObject    handle to Load_data_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
[fnames,fpath] = uigetfile('*Ch1.txt','multiselect','on','Select PDV files');
set(handles.text_prompt,'String','Loading, Please Wait');
if fpath ==0
    error('LoadData:NoData',...
        'No file input')
end

if ischar(fnames)
    fnames = {fnames};
end
[handles.time,handles.amplitude,handles.t0,handles.time_offset] = file_read(fpath,fnames,handles.scope_offset);
[handles.camp] = detrend_amp(hObject,eventdata,handles);
%At this point we cand delete old data
Delete_Old_Data(hObject,eventdata,handles);
handles.fpath = fpath; handles.fnames = fnames; handles.file_count = length(fnames);
set(handles.files_listbox,'Value',1); 
set(handles.files_listbox,'String',[handles.fnames']);
handles.list_idx = handles.file_count;
handles.current_message = '';
handles.calc_bool = 0; handles.peaks_bool = 0;
handles.rise_t_bool = 0; set(handles.rise_t_button,'Value',0);
guidata(hObject,handles);
plot_int(hObject,eventdata,handles);
error('Allfunctions:NoError','No Error');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end
% --- Executes on selection change in files_listbox.
function files_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns files_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from files_listbox
try
    if isempty(handles.fnames)
        error('LISTBOX:NOFILE','No Files Loaded');
    end     
handles.list_idx = get(hObject,'Value');
guidata(hObject,handles);
switch handles.rise_t_bool
    case 0
        update_plots(hObject,eventdata,handles);
    case 1
        update_plots_IM(hObject,eventdata,handles);
end
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end


% --- Executes during object creation, after setting all properties.
function files_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_all_button.
function run_all_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_all_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
if isempty(handles.fnames)
    error('RunData:NoData','No Data to Run');
end

for i = 1:handles.file_count
    handles.list_idx = i;
    handles.current_message = sprintf('Processed %s',...
        handles.fnames{handles.list_idx});
    switch handles.method
        case 1
             handles.phase{i} = phase_analysis(handles.camp{i});
            [handles.lineout_time{i},handles.velocity{i},handles.min_idx{i},handles.displacement{i}] = STFT_analysis(hObject,eventdata,handles);
            handles.velocity{i} = handles.velocity{i}./handles.window_val;
            [handles.pressure{i},handles.flux{i},handles.fluence{i}] = calc_derived_data(hObject,eventdata,handles);
            handles.peaks_bool =0;
            update_plots(hObject,eventdata,handles);
        case 2
            handles.phase{i} = phase_analysis(handles.camp{i});
            [~,handles.min_idx{i}] = min(handles.phase{i});
            [handles.lineout_time{i},handles.velocity{i},handles.displacement{i},handles.xyPeaks{i}] = peak_det_3(hObject,eventdata,handles);
            handles.velocity{i} = handles.velocity{i}./handles.window_val;
            handles.peaks_bool =1;
            [handles.pressure{i},handles.flux{i},handles.fluence{i}] = calc_derived_data(hObject,eventdata,handles);
            %handles.xPeaks{i} = handles.xypeaks{i}(:,1);handles.yPeaks{i} = handles.xypeaks{i}(:,2);
            update_plots(hObject,eventdata,handles);
    end
end
handles.calc_bool =1;
guidata(hObject,handles);
handles.current_message = '';
guidata(hObject,handles);
error('Allfunctions:NoError','No Errors');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end

% --- Executes on button press in phase_min_button.
function phase_min_button_Callback(hObject, eventdata, handles)
% hObject    handle to phase_min_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.lineout_axes);
n = handles.list_idx;
[x,~] = ginput(1);
switch handles.rise_t_bool
    case 0
        [~,handles.min_idx{n}] = min(abs(handles.time{n}-x));
    case 1
        [~,handles.min_idx{n}] = min(abs((handles.time{n}-handles.time{n}(handles.t0{n}))-x));
end
[handles.velocity{n}] = update_phase_min(hObject,eventdata,handles);
[handles.pressure{n},handles.flux{n},handles.fluence{n}] = calc_derived_data(hObject,eventdata,handles);
[handles.displacement{n}] = update_displacement(hObject,eventdata,handles);
switch handles.rise_t_bool
    case 0
        update_plots(hObject,eventdata,handles);
    case 1
        update_plots_IM(hObject,eventdata,handles);
end
guidata(hObject,handles);


% --- Executes on selection change in plot_type_list.
function plot_type_list_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_type_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type_list
handles.plot_idx = get(hObject,'Value');
update_plots(hObject,eventdata,handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function plot_type_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_tab_Callback(hObject, eventdata, handles)
% hObject    handle to File_tab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function settings_menu_Callback(hObject, eventdata, handles)
% hObject    handle to settings_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    PDV_flux_3_Settings;
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end

% --------------------------------------------------------------------
function add_peaks_tab_Callback(hObject, eventdata, handles)
% hObject    handle to add_peaks_tab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    PDV_flux_3_Change_Peaks; %uiwait(PDV_flux_3_Change_Peaks);
    error('Allfunctions:NoError','No Error');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end
    

% --------------------------------------------------------------------
function Save_data_item_Callback(hObject, eventdata, handles)
% hObject    handle to Save_data_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if handles.calc_bool ==0
        error('SaveFunction:NoSaveData','Nothing to save')
    end
hdr1 = {}; hdr2 = {}; hdr3 = {};
max_vector_size = [];
    for i=1:handles.file_count
        max_vector_size(i) = length(handles.lineout_time{i});
    end
    max_vector_size = max(max_vector_size);
full_save = {};
for i = 1:handles.file_count
    curr_size = length(handles.lineout_time{i});
    save_data = [];
    save_data(1:max_vector_size,1:2) = NaN;
    hdr1 = [hdr1,'Time','Velocity'];
    hdr2 = [hdr2,'ns','km/s'];
    hdr3 = [hdr3, handles.fnames{i},handles.fnames{i}];
    save_data(1:curr_size,:) = [handles.lineout_time{i},handles.velocity{i}];
    full_save{i} = save_data;
end
if ~isempty(full_save)
    work_dir = pwd;
    cd(handles.fpath); filter = {'*.txt'};
    [save,save_path] = uiputfile(filter,'Save PDV file');
    if save_path ==0
        cd(work_dir);
        error('SaveFunc:CancelInput','Save Cancelled');
    end
    fmt = repmat('%s\t ', 1, length(hdr1));
    fmt(end:end+1) = '\n';
    %open save file and write headers
    fid = fopen(fullfile(save_path,save), 'w');
    fprintf(fid, fmt, hdr1{:});
    fprintf(fid,fmt, hdr2{:});
    fprintf(fid,fmt, hdr3{:});
    fclose(fid);
    %now insert data vector
    dlmwrite(fullfile(save_path,save),full_save,'-append','delimiter','\t');
    handles.current_message = 'File Saved';
    cd(work_dir);
else
    handles.current_message = 'No files to save';
end
error('Allfunctions:NoError','No Error');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end

function TIM_menu_Callback(hObject, eventdata, handles)
% hObject    handle to TIM_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PDV_flux_3_TIM;

% --- Executes on button press in rise_t_button.
function rise_t_button_Callback(hObject, eventdata, handles)
% hObject    handle to rise_t_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.rise_t_bool = get(hObject,'Value');
switch handles.rise_t_bool
    case 1
        update_plots_IM(hObject,eventdata,handles);
    case 0
        update_plots(hObject,eventdata,handles);
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function pressure_save_Callback(hObject, eventdata, handles)
% hObject    handle to pressure_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles.save_type = 1;
save_der_data(hObject,eventdata,handles)
guidata(hObject,handles);
handles.current_message = 'Data Saved';
error('Allfunctions:NoError','No Error');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end

% --------------------------------------------------------------------
function flux_save_Callback(hObject, eventdata, handles)
% hObject    handle to flux_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles.save_type = 2;
save_der_data(hObject,eventdata,handles)
guidata(hObject,handles);
handles.current_message = 'Data Saved';
error('Allfunctions:NoError','No Error');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end

% --------------------------------------------------------------------
function disp_save_Callback(hObject, eventdata, handles)
% hObject    handle to disp_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles.save_type = 3;
save_der_data(hObject,eventdata,handles)
guidata(hObject,handles);
handles.current_message = 'Data Saved';
error('Allfunctions:NoError','No Error');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end

% --------------------------------------------------------------------
function fluence_save_Callback(hObject, eventdata, handles)
% hObject    handle to fluence_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles.save_type = 4;
save_der_data(hObject,eventdata,handles)
guidata(hObject,handles);
handles.current_message = 'Data Saved';
error('Allfunctions:NoError','No Error');
catch ME
    exception_handler(ME,hObject,eventdata,handles);
end

% --------------------------------------------------------------------
function save_options_Callback(hObject, eventdata, handles)
% hObject    handle to save_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to flux_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function TIM_save_Callback(hObject, eventdata, handles)
% hObject    handle to TIM_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.launch_save,'Checked','off');
set(hObject,'Checked','on');
handles.save_timing_bool = 1; %means it will plot starting at rise
guidata(hObject,handles);
% --------------------------------------------------------------------
function launch_save_Callback(hObject, eventdata, handles)
% hObject    handle to launch_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TIM_save,'Checked','off');
set(hObject,'Checked','on');
handles.save_timing_bool = 0; %means it will plot starting at rise
guidata(hObject,handles);

% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to disp_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------


% --------------------------------------------------------------------
function save_settings_Callback(hObject, eventdata, handles)
% hObject    handle to save_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_timing_select_Callback(hObject, eventdata, handles)
% hObject    handle to save_timing_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in reset_t0_button.
function reset_t0_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_t0_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.lineout_axes);
n = handles.list_idx;
[x,~] = ginput(1);
switch handles.rise_t_bool
    case 0
        [~,new_t0] = min(abs(handles.time{n}-x));
    case 1
        [~,new_t0] = min(abs((handles.time{n}-handles.time{n}(handles.t0{n}))-x));
end
handles.t0{n} = new_t0;
[handles.pressure{n},handles.flux{n},handles.fluence{n}] = calc_derived_data(hObject,eventdata,handles);
[handles.displacement{n}] = update_displacement(hObject,eventdata,handles);
switch handles.rise_t_bool
    case 0
        update_plots(hObject,eventdata,handles);
    case 1
        update_plots_IM(hObject,eventdata,handles);
end
guidata(hObject,handles);
%%Begin subroutine section

function [time,amplitude,t0,time_offset] = file_read(fpath,fnames,scope_offset)
    [time,amplitude,t0,time_offset] = cellfun(@(x) channel_read(fpath,x,scope_offset),fnames,'UniformOutput',false);

%% reads enumerated list of files
function [time,amplitude,t0,time_offset] = channel_read(fpath,fname,scope_offset)
    name = strsplit(fname,'Ch'); name = name{1};
    for i = 1:4
        channel_name=strcat(name,sprintf('Ch%d.txt',i));
        fid = fopen(fullfile(fpath,channel_name));
        textscan(fid,'%s',31);
        file = fscanf(fid,'%f',[2,1])';
        while ~feof(fid)
            curr = fscanf(fid,'%f',[2,5000])';
            if ~isempty(curr)
                file = [file; curr];
            end
        end 
        fclose(fid);
        if i==1
            time = file(:,1) + abs(file(1,1));
        end
        amplitude(:,i) = file(:,2);
        clear file
        s=0;
        for j=1:(0.05*length(time))
            s=s+amplitude(j,i)^2;
        end
        rMS(i)=sqrt(s/(0.05*length(time)));

    end
    
    i = 1;
    try
    while (abs(amplitude(i,1)) < 5*rMS(1) || abs(amplitude(i,2)) < 5*rMS(2) && i < length(amplitude(:,1)))
        i = i+1;
    end
    catch
        i = 1000;
    end
    [maximum, maximum_index] = max(amplitude(:,3));
    time_vector = time(1:maximum_index);
    index90 = length(time_vector(time_vector<=maximum*0.9));
    time90 = time(index90).*1e9;
    time_offset = -time90 + scope_offset;
    time = time.*1e9 + time_offset;
    t0 = i;
%% Detrending of amplitude files and normalizing
    function [camp] = detrend_amp(hObject,eventdata,handles)
        [camp] = cellfun(@(x,y,z) corr_amp(x,y,z), handles.time,handles.amplitude,handles.t0,'UniformOutput',false);
        
    function [camp] = corr_amp(time,amp,t0) 
        A = [0.35,0.37,1,0.34];
        for i = 1:4
            if i ==3
                camp(:,i) = amp(:,i);
                continue
            end
            fitvals = polyfit(time(t0:end),amp(t0:end,i),1);
            meanvals = polyval(fitvals,time);
            camp(:,i) = (amp(:,i)-meanvals)./A(i);
        end

%% central finite difference derivative (first order)
function [der] = dx_dt(time,amplitude)
    for k=1:length(amplitude(1,:))
        der(1,k) = 0;
        for i=2:length(time)-1
            der(i,k) = (amplitude(i+1,k)-amplitude(i-1,k))/(time(i+1)-time(i-1)); 
        end
        der(i+1,k) = 0;
    end
%% central finite difference derivative (second order)
function [der2] = dx_dt_2(time,amplitude)
    for k=1:length(amplitude(1,:))
        amplitude(:,k) = smooth(amplitude(:,k),5);
        dt = time(2)-time(1);
        der2(1:2,k) = 0;
        for i=3:length(time)-2
            der2(i,k) = (-amplitude(i+2,k)+8*amplitude(i+1,k)-8*amplitude(i-1,k)+amplitude(i-2,k))/(12*dt);
        end
        der2(i+1:i+2,k) = 0;
    end

%% peak finding function, computes derivative then finds zero crossings
    function [lineout_time,velocity,xpeaks,ypeaks] = Fringe_analysis(hObject,eventdata,handles)
            n = handles.list_idx;
            [xpeaks,idx_peaks,lineout_time,velocity] = find_peaks(handles.time{n},handles.amplitude{n},handles.t0{n});
            amp = handles.amplitude{n}; 
            for i=1:4; ypeaks{i} = amp(idx_peaks{i},i); 
            end
            
        
        
    function [xpeaks,idx_peaks,lineout_time,velocity] = find_peaks(time,amplitude,t0)
    der = dx_dt(time,amplitude);
    [xpeaks,idx_peaks] = find_zeros(time,der,t0);
    xpeaks{1} = sortrows(xpeaks{1});xpeaks{2} = sortrows(xpeaks{2});xpeaks{4} = sortrows(xpeaks{4});
    velocity1 = [0;0.3875./diff(xpeaks{1})];velocity2 = [0;0.3875./diff(xpeaks{2})];
    velocity3 = [0;0.3875./diff(xpeaks{4})];
    lineout_time = sortrows([xpeaks{1},velocity1;xpeaks{2},velocity2;xpeaks{4},velocity3]);
    velocity = lineout_time(:,2); lineout_time = lineout_time(:,1);
    
    
%% finds the zero crossing of the periodic data function. finds first where the sign changes then draws a
%order 1 polynomial between and finds the root
function [xpeak,ypeak] = find_zeros(x,y,x0)

    for j=1:4
        %y(:,j) =smooth(y(:,j),'sgolay',3);
        %y(:,j) = smooth(y(:,j),3);
        locs = [];
        k = x0-4;
        if isempty(k)
            error('Bad x0 value')
        end
        init  = 0;
        while init ==0
            if y(k,j)<0
                locs = [locs,x(k)];
                init = -1;
            elseif y(k,j)>0
                locs = [locs,x(k)];
                init = 1;
            else
                k = k+1;
            end
        end
        i = k;
        while i<length(x)
            while init>0 && i<length(x)
                %look for first negative, then linearly interpolate the zero
                %and index it. We assume the the velocity is below ~4 so it
                %imposes that the time spacing must be >0.1 ns and <100 ns
                if y(i,j)<0
                    p = polyfit(x(i-1:i),y(i-1:i,j),1);
                    zer_loc = roots(p);
                    if abs(locs(end)-zer_loc)>0.2 && abs(locs(end)-zer_loc)<100
                            locs = [locs,zer_loc];
                            init = -1;
                            
                    else
                        i=i+1;
                    end
                elseif y(i,j)>=0
                    i = i+1;
                    
                end
            end
            
                while init <0 && i<length(x)
                    %look for first positive, then linearly interpolate the
                    %zero
                    if y(i,j)>0
                        p = polyfit(x(i-1:i),y(i-1:i,j),1);
                        zer_loc = roots(p);
                        if abs(locs(end)-zer_loc)>0.2 && abs(locs(end)-zer_loc)<100
                            
                                locs = [locs,zer_loc];
                                init = 1;
                                
                            
                        else
                            i=i+1;
                        end
                    elseif y(i,j)<=0
                        i = i+1;
                        
                    end
                end
            
        end
        xpeaks{j} = locs';
        clear ypeaks
        for k=1:length(locs)
            idx = find(abs(x-locs(k))==min(abs(x-locs(k))));
            ypeaks(k) = idx(1);
        end
        ypeak{j} = ypeaks';
    end
    xpeak = {xpeaks{1},xpeaks{2},xpeaks{3},xpeaks{4}};
    ypeak = {ypeak{1},ypeak{2},ypeak{3},ypeak{4}};
%% Short time fourier transform analysis

    function [lineout_time,velocity_lineout_fit,minimum_index,displacement] = STFT_analysis(hObject,eventdata,handles)
        n = handles.list_idx;
        time = handles.time{handles.list_idx};
        sample_frequency = (time(2)-time(1))*1e-9;
        sample_frequency = sample_frequency^(-1);
        sample_spacing = time(2)-time(1);
        r = round(sample_frequency.*(handles.time_resolution*1e-9));
        test = round(handles.sample_rate/sample_spacing);
        camp = handles.camp{handles.list_idx};
            if test==0
                test=1;
            end
        for i = 1:4
            if i==3
                continue
            end
            [STFT,f,t] = spectrogram(camp(:,i),hamming(r),r-test,10*r,sample_frequency);
            if i==1
                STFT_tot = abs(STFT);
            else
                STFT_tot = STFT_tot + abs(STFT);
            end
        end
        STFT_tot = STFT_tot./3;
           
            velocity_axis = f.*0.775./1e9;
            lineout_time = (t*1e9 + handles.time_offset{handles.list_idx})';
            if handles.velocity_cutoff_low > 0
                filter = length(velocity_axis(velocity_axis<handles.velocity_cutoff_low));
                STFT_tot(1:filter,:)=0;
                clear filter
            end
            [mx locs]=max(STFT_tot,[],1);
            velocity_lineout=velocity_axis(locs);

            % Fit the FFT at each time step to better resolve the velocity. I use a
            % polynomial since this is much, much less computationally expensive then a
            % gaussian fit. 
            velocity_lineout_fit = velocity_lineout;
            for i=1:length(velocity_lineout)
                if velocity_lineout(i) > 0.1 && (locs(i)+2)<length(velocity_axis)
                    p = polyfit(velocity_axis((locs(i)-2):(locs(i)+2)),STFT_tot((locs(i)-2):(locs(i)+2),i),2);
                    peakPosition = -p(2)./(p(1)*2);
                    velocity_lineout_fit(i) = peakPosition;
                else
                    velocity_lineout_fit(i) = 0;
                end
            end
            [~, minimum_index] = min(handles.phase{handles.list_idx});
            displacement = [];
for i = 1:length(lineout_time)
    if lineout_time(i) < handles.time{n}(handles.t0{n})
        displacement(i) = 0;
    else
        displacement(i) = displacement(i-1)+trapz(lineout_time(i-1:i),velocity_lineout_fit(i-1:i));
    end
end
for i=1:length(lineout_time)
    if i==1
        velocity_lineout_fit(i) = velocity_lineout_fit(i);
    elseif lineout_time(i)<time(minimum_index)
        velocity_lineout_fit(i) = velocity_lineout_fit(i);
    else
        velocity_lineout_fit(i) = -1*velocity_lineout_fit(i);
    end
end


        
guidata(hObject,handles);
            
         
            
%% Old peak finding program, credit to Will Shaw and his undergrad.            
function [lineout_time,velocity_final,displacement_s,xyPeaks] = peak_det_3(hObject,eventdata,handles)
 n = handles.list_idx;
 amp = handles.amplitude{n};
 time = handles.time{n};
 x0 = handles.t0{n}-4;
 time0 = time(x0);
 amp0 = amp(x0,1);
 k = 1;
 for i = 1:4
     if i==3
         continue
     end
     samp(:,k) = smooth(amp(:,i),5);
     rMS(k) = fRMS(time,amp(:,i));
     [high{k},low{k}] = peakdet(samp(x0-0.001*length(time):end,k),(rMS(k)*handles.threshold_val),time(x0-0.001*length(time):end));
     peakPositions{k} = sortrows([time(x0),amp(x0,i);high{k}(:,1),high{k}(:,2);low{k}(:,1),low{k}(:,2)]);
     [samp_test] = fSmoothData2(time,amp(:,i),peakPositions{k},0.1);
     if length(samp_test) == length(time)
        samp(:,k) = samp_test;
     end
     [high{k},low{k}] = peakdet(samp(x0:end,k),(rMS(k)*handles.threshold_val),time(x0:end));
     peakPositions{k} = sortrows([time(x0),amp(x0,i);high{k}(:,1),high{k}(:,2);low{k}(:,1),low{k}(:,2)]);
     %Call the fitting program that uses a third 2nd order polynomial to fit the
     %data about each peak found above. This gives the most precise position for
     %each max or min.
     [high{k},low{k}] = fFitPeaks(time,amp(:,i),high{k},low{k});
     peakPositions{k} = sortrows([time(x0),amp(x0,i);high{k}(:,1),high{k}(:,2);low{k}(:,1),low{k}(:,2)]);
     xyPeaks{k}=sortrows([time(x0),samp(x0,k);high{k}(:,1),high{k}(:,2);low{k}(:,1),low{k}(:,2)]);
     idx_0 = length(xyPeaks{k}(xyPeaks{k}(:,1)<=time(x0),1));
     velocity{k}=0.3875./diff(xyPeaks{k}(idx_0+1:end,1));
     xPeaks{k}=xyPeaks{k}(idx_0+1:end,1);xPeaks{k}(length(velocity{k}))=[];
     k = k+1;
 end 
    velocity0 = 0;
    XYMAT=sortrows([xPeaks{1},velocity{1};xPeaks{2},velocity{2};xPeaks{3},velocity{3};time(x0),velocity0]);
    velocityTime=XYMAT(:,1);
    velocity_final=XYMAT(:,2);
    x = smooth(velocity_final(3:end),3);
    velocity_final(3:end) = x;
    lineout_time = XYMAT(:,1);
    new_time = time-velocityTime(1,1);
for k = 1:3
peakPositions1{k}(:,1) = peakPositions{k}(:,1)-velocityTime(1,1);
xPeaks{k} = xPeaks{k} - velocityTime(1,1);
end

displacement1=cell(length(xPeaks{1}),1);
displacement2=cell(length(xPeaks{2}),1);
displacement3=cell(length(xPeaks{3}),1);
[~, minimum_index] = min(handles.phase{n});

for i=1:length(lineout_time)
    if i==1
        velocity_final(i) = velocity_final(i);
    elseif lineout_time(i)<time(minimum_index)
        velocity_final(i) = velocity_final(i);
    else
        velocity_final(i) = -1*velocity_final(i);
    end
end
for i=1:length(xPeaks{1})
    if i==1
        displacement1{i}=(1.55/(4*1.0627));
    elseif xPeaks{1}(i) < time(minimum_index)
        displacement1{i}=displacement1{i-1}+(1.55/(4*1.0627));
    else
        displacement1{i}=displacement1{i-1}-(1.55/(4*1.0627));
    end
end
for i=1:length(xPeaks{2})
    if i==1
        displacement2{i}=(1.55/(4*handles.window_val));
    elseif xPeaks{2}(i) < time(minimum_index)
        displacement2{i}=displacement2{i-1}+(1.55/(4*handles.window_val));
    else
        displacement2{i}=displacement2{i-1}-(1.55/(4*handles.window_val));
    end
end
for i=1:length(xPeaks{3})
    if i==1
        displacement3{i}=(1.55/(4*1.0627));
    elseif xPeaks{3}(i) < time(minimum_index)
        displacement3{i}=displacement3{i-1}+(1.55/(4*handles.window_val));
    else
        displacement3{i}=displacement3{i-1}-(1.55/(4*handles.window_val));
    end
end
displacement_c=cell(3,2);
displacement_c{1,1}=xPeaks{1};displacement_c{2,1}=xPeaks{2};displacement_c{3,1}=xPeaks{3};
displacement_c{1,2}=cell2mat(displacement1);displacement_c{2,2}=cell2mat(displacement2);displacement_c{3,2}=cell2mat(displacement3);
displacement2c=cell2mat(displacement_c);
displacement2c=sortrows(displacement2c);
displacementTime = displacement2c(:,1);
displacement=displacement2c(:,2);
displacement_s=smooth(displacementTime,displacement,0.1,'rloess');
displacement_s = [0;displacement_s];
guidata(hObject,handles);
%Calculate final time to end integration
% j=2;
% while j<length(flux)&&(fluxTime(j,1)-fluxTime((j-1),1))<5
%     j=j+1;
% end
%}

%% Sub Functions for old peakdet program
function [rMS] = fRMS(time,rVolts1)
s=0;
for i=1:(0.02*length(time))
    s=s+rVolts1(i,1)^2;
    i=i+1;
end
rMS=sqrt(s/(0.02*length(time)));


%%
function [sVolts] = fSmoothData2(rTime,rVolts1,peakPositions,num)
peakPositionsIndex = [];
for i=1:length(peakPositions(:,1))
    peakIndex = find(rTime >= peakPositions(i,1),1,'first');
    peakPositionsIndex = [peakPositionsIndex;peakIndex];
end

indexDifference = [peakPositionsIndex(1);diff(peakPositionsIndex)];


sVolts = smooth(rVolts1(1:peakPositionsIndex(1)),(indexDifference(1)*num));

smoothDistance = round((indexDifference(2)*num));
z = 2;

while smoothDistance<1
    smoothDistance = round((indexDifference(z+1)*num));
    z = z+1;
end

for i = 2:(length(peakPositions(:,1))-10)
    if round((indexDifference(i)*num)) > (1+1/4)*smoothDistance
        smoothDistance = round((1+1/4)*smoothDistance);
    elseif round((indexDifference(i)*num)) < (1-1/4)*smoothDistance
        smoothDistance = round((1-1/4)*smoothDistance);
    else
        smoothDistance = round((indexDifference(i)*num));
    end
    if i==length(peakPositions(:,1))
        disp('here')
    end
    test = smooth(rVolts1((peakPositionsIndex(i-1)-smoothDistance):(peakPositionsIndex(i)+smoothDistance)),smoothDistance);
    sVolts = [sVolts;test(smoothDistance+1:end-smoothDistance-1)];
end
test = smooth(rVolts1((peakPositionsIndex(end)+1)-100:end),(50));
sVolts = [sVolts;test(101:end)];


%%
function [newHigh,newLow] = fFitPeaks(xData,yData,high,low)
%Title: findPeaks
%Author: Gino Giannetti & William Shaw

%Date: 2014-6-10
%Updated 2014-12-31

%Purpose: find maximums and minimums of input data more accurately than
%         peakdet method using curve fitting. For use in shock compression data.
%         Note-uses peakdet in the process (Removed dependency 2014-12-31)

%Input: high and low are previous max and min arrays found by peakdet. 
%       xData and yData are the raw x and y data 
%       indFirstMove is the index in the x and y data where the object in
%       question first moves

%Output: new maximum and minimum arrays of the data

% Setting up while loop

oldPeaks = sortrows([high(:,1),high(:,2);low(:,1),low(:,2)]);
newHigh = [];
newLow = [];
a=0;
b=0;
k=1;

warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
while (k <= (length(oldPeaks(:,1))))

% Finding section to curve fit
    switch k
        case 1
            a = ((oldPeaks(k,1)-oldPeaks(k+1,1)))/2 + oldPeaks(k,1);
            b = ((oldPeaks(k+1,1)-oldPeaks(k,1)))/2 + oldPeaks(k,1); 
        case length(oldPeaks(:,1))
            a = ((oldPeaks(k,1)-oldPeaks(k-1,1)))/2 + oldPeaks(k-1,1);
            b = ((oldPeaks(k,1)-oldPeaks(k-1,1)))/2 + oldPeaks(k,1);
        otherwise
            a = ((oldPeaks(k,1)-oldPeaks(k-1,1)))/2 + oldPeaks(k-1,1);
            b = ((oldPeaks(k+1,1)-oldPeaks(k,1)))/2 + oldPeaks(k,1);
            
    end
    
    c=length(xData(xData<=a));
    d=length(xData(xData<=b));
    
    if c<0
        c=1;
    end
    
    if(isempty(c) || isempty(d))
        display('peak could not be found');
    end
    
    if c > length(xData) || d > length(xData)
        k=k+1;
    else
        xTemp = xData((c):(d));
        yTemp = yData((c):(d));
    
        % curve fitting polynomial 
        p = polyfit(xTemp,yTemp,2);
        peakPosition = -p(2)./(p(1)*2);
        if peakPosition > xData(end)
            peakIndex = 1;
        else
            peakIndex = find(xData >= peakPosition,1,'first');
        end
        
        % Don't allow peak positions to change by more than 1 ns
        positionChange = peakPosition-oldPeaks(k,1);
        if abs(positionChange) < 1
            Temp = [peakPosition,yData(peakIndex)];
        else
            Temp = [oldPeaks(k,1), oldPeaks(k,2)];
        end
        
        %Assign high or low value.
        if (p(1)<=0)
            newHigh = [newHigh;Temp];
        else 
        newLow = [newLow;Temp];
        end
    
        k=k+1;
    end
end
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');

%% Peak determining routine
function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.

maxtab = [];
mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
  x = (1:length(v))';
else 
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end
  
if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(v)
  this = v(i);
  if this > mx, mx = this; mxpos = x(i); end
  if this < mn, mn = this; mnpos = x(i); end
  
  if lookformax
    if this < mx-delta
      maxtab = [maxtab ; mxpos mx];
      mn = this; mnpos = x(i);
      lookformax = 0;
    end  
  else
    if this > mn+delta
      mintab = [mintab ; mnpos mn];
      mx = this; mxpos = x(i);
      lookformax = 1;
    end
  end
end

%% Complex analysis functions
function [phase] = phase_analysis(camp)
%Correct this later for 3-channel.
%Form the complex interferrogram
mixX = sin(-120*(pi/180))*camp(:,1);                %Correct Equation
mixY = cos(-120*(pi/180))*camp(:,1)-camp(:,2);          %Correct Equation
complexAmp = mixX - 1i.*mixY;

%Calulate the phase angle change. The minimum corresponds to the change in
%direction.
phase = unwrap(angle(complexAmp));
%[minimum, minimum_index] = min(phase);


%Find the absolute time offset based on trigger 3
%[z,zloc]=max(rVolts3); %identify the maximum value
%{
Leaving this out until I figure out how I want complex spectra to be
deployed i UI
sden2=rVolts3(1:zloc); %remove all data after maximum value
hm=length(sden2(sden2<=z.*0.9));  %identify value of 90% for time 0
z=time(hm);   %report correction for time 0 value in nanoseconds
tAbsolute=-z+12;    %time correction offset
%}
%% Function to calculate various kinds of derived data

    function [pressure,flux,fluence] = calc_derived_data(hObject,eventdata,handles)
        i=1;
        j=1;
        n = handles.list_idx;
        time = handles.time{n};
        velocity_lineout_fit = handles.velocity{n};
        time_lineout = handles.lineout_time{n};
        fluence(j,1) = 0;
        pressure(j,1)=0;
        x = 1;
        [~,minimum_index] = min(handles.phase{n});
        switch handles.sample_ID
            case 1
            rho = 2.230;
        while i<length(velocity_lineout_fit)
            j=j+1;
            if time_lineout(i)> time(minimum_index)
                        x = -1;
            end
            if abs(velocity_lineout_fit(i,1)) <= 0.568
                m = 1.861;%% note: these are specific to GLASS, from the glass hugoniot --> that is why the constants change for different velocities 
                b = 3.879;
                flux(j,1)=find_flux(m,b,rho,velocity_lineout_fit(i,1));% the flux is = force * dA , where A is area 
                pressure((j),1)=find_pressure(m,b,rho,velocity_lineout_fit((i),1));
                fluence(j,1)=fluence((j-1),1)+x*trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
            elseif abs(velocity_lineout_fit(i,1)) > 1.83
                m = 1.269;
                b = 2.3925;
                flux(j,1)=find_flux(m,b,rho,velocity_lineout_fit(i,1)); 
                pressure((j),1)=find_pressure(m,b,rho,velocity_lineout_fit((i),1));
                fluence(j,1)=fluence((j-1),1)+x*trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
            else 
                m = -0.175;
                b = 5.0344; 
                flux(j,1)=find_flux(m,b,rho,velocity_lineout_fit(i,1));% the flux is = force * dA , where A is area 
                pressure((j),1)=find_pressure(m,b,rho,velocity_lineout_fit((i),1));
                fluence(j,1)=fluence((j-1),1)+x*trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
            end 
            i=i+1;
        end
            case 2 %Nitromethane
                rho = 1.125; m = 1.64; b = 1.65;
                while i<length(velocity_lineout_fit)
                    j=j+1;
                    if time_lineout(i)> time(minimum_index)
                        x = -1;
                    end
                    flux(j,1) = find_flux(m,b,rho,velocity_lineout_fit(i,1));
                    pressure(j,1) = find_pressure(m,b,rho,velocity_lineout_fit(i,1));
                    fluence(j,1) = fluence((j-1),1)+x*trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
                    i = i+1;
                end
            case 3 %Sapphire
                rho = 4.0; m = 0.957; b = 8.74;
                while i<length(velocity_lineout_fit)
                    j=j+1;
                    if time_lineout(i)> time(minimum_index)
                        x = -1;
                    end
                    flux(j,1) = find_flux(m,b,rho,velocity_lineout_fit(i,1));
                    pressure(j,1) = find_pressure(m,b,rho,velocity_lineout_fit(i,1));
                    fluence(j,1) = fluence((j-1),1)+x*trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
                    i = i+1;
                end
            case 4 %CaF2
                %Double check this
                rho = 3.18; m = 1.18; b = 5.15;
                while i<length(velocity_lineout_fit)
                    j=j+1;
                    if time_lineout(i)> time(minimum_index)
                        x = -1;
                    end
                    flux(j,1) = find_flux(m,b,rho,velocity_lineout_fit(i,1));
                    pressure(j,1) = find_pressure(m,b,rho,velocity_lineout_fit(i,1));
                    fluence(j,1) = fluence((j-1),1)+x*trapz(time_lineout((j-1):j,1),flux((j-1):j,1));
                    i = i+1;
                end
        end
        
%% Find flux and fluence equations;        
    function flux = find_flux(m,b,rho,velocity)
        flux = 0.5*rho*(b+m*velocity)*velocity^2; % the flux is = force * dA , where A is area


%%function to calculate pressure given m,p,b,and vel
    function pressure = find_pressure(m,b,rho,velocity)
        pressure = rho*(b+m*velocity)*velocity;
%% Saving Derived Data

    function save_der_data(hObject,eventdata,handles)
    if handles.calc_bool ==0
        error('SaveFunction:NoSaveData','Nothing to save')
    end
    hdr1 = {}; hdr2 = {}; hdr3 = {};
    max_vector_size = [];
    for i=1:handles.file_count
        max_vector_size(i) = length(handles.lineout_time{i});
    end
    max_vector_size = max(max_vector_size);
    full_save = {};
    switch handles.save_type
        case 1 %Pressure
            x = 'pressure'; y = 'GPa';
            der_data = handles.pressure;
        case 2 %Flux
            x = 'flux'; y = 'TJ/m-s';
            der_data = handles.flux;
        case 3 %Displacement
            x = 'displacement'; y = 'um';
            der_data = handles.displacement;
        case 4 %Fluence
            x = 'fluence'; y = 'kJ/m^2';
            der_data = handles.fluence;
    end
    for i = 1:handles.file_count
        switch handles.save_timing_bool
            case 0
                time = handles.lineout_time{i};
            case 1
                time = handles.lineout_time{i}-handles.time{i}(handles.t0{i});
        end
        curr_size = length(handles.lineout_time{i});
        save_data = [];
        save_data(1:max_vector_size,1:3) = NaN;
        hdr1 = [hdr1,'time','velocity',x];
        hdr2 = [hdr2,'ns','km/s',y];
        hdr3 = [hdr3, handles.fnames{i},handles.fnames{i},handles.fnames{i}];
        save_data(1:curr_size,:) = [time,handles.velocity{i},der_data{i}];
        full_save{i} = save_data;
    end
    if ~isempty(full_save)
        work_dir = pwd;
        cd(handles.fpath); filter = {'*.txt'};
        [save,save_path] = uiputfile(filter,'Save PDV file');
        if save_path ==0
            cd(work_dir);
            error('SaveFunc:CancelInput','Save Cancelled');
        end
        fmt = repmat('%s\t ', 1, length(hdr1));
        fmt(end:end+1) = '\n';
        %open save file and write headers
        fid = fopen(fullfile(save_path,save), 'w');
        fprintf(fid, fmt, hdr1{:});
        fprintf(fid,fmt, hdr2{:});
        fprintf(fid,fmt, hdr3{:});
        fclose(fid);
        %now insert data vector
        dlmwrite(fullfile(save_path,save),full_save,'-append','delimiter','\t');
        handles.current_message = 'File Saved';
        cd(work_dir);
    else
        error('SaveFunction:NoFiles', 'No files to save');
    end
%%updating phase minimum
        function new_velocity = update_phase_min(hObject,eventdata,figure1)
            n = figure1.list_idx;
            figure1.velocity{n} = abs(figure1.velocity{n});
            for i = 1:length(figure1.velocity{n})
                if figure1.lineout_time{n}(i)<figure1.time{n}(figure1.min_idx{n})
                    figure1.velocity{n}(i) = figure1.velocity{n}(i);
                else
                    figure1.velocity{n}(i) = -1.*figure1.velocity{n}(i);
                end
            end
            new_velocity = figure1.velocity{n};
            
            function new_displacement = update_displacement(hObject,eventdata,handles)
                n = handles.list_idx;
                vel = handles.velocity{n};
                l_time = handles.lineout_time{n};
               
                switch handles.peaks_bool
                    case 0  %STFT USED
                         new_displacement = [];
                         for i = 1:length(l_time)
                            if l_time(i) < handles.time{n}(handles.t0{n})
                                new_displacement(i) = 0;
                            else
                                new_displacement(i) = new_displacement(i-1)+trapz(l_time(i-1:i),vel(i-1:i));
                            end
                        end
                    case 1 %Peak finding used
                        [~,~,new_displacement,~] = peak_det_3(hObject,eventdata,handles);
                end
            
%% plots interferogram
    function int = plot_int(hObject,eventdata,handles)
        n=handles.list_idx;
        amp = handles.amplitude{n};
        axes(handles.interferogram_axes); hold off;
        for i = 1:4
            if i==3
                continue
            end
            plot(handles.time{n},amp(:,i)); hold on;
            if (handles.method ==2 && handles.peaks_bool ==1)
                if i==4
                    xvals = handles.xyPeaks{n}{3}(:,1);yvals = handles.xyPeaks{n}{3}(:,2);
                else
                    xvals = handles.xyPeaks{n}{i}(:,1);yvals = handles.xyPeaks{n}{i}(:,2);
                end
                plot(xvals,yvals,'ok')
            end
        end
        xlabel('time (ns)');ylabel('volts (V)');
        if handles.peaks_bool ==0
            legend('Channel 1','Channel 2','Channel 3');
        else
            legend('Channel 1','','Channel 2','','Channel 3','');
        end
        
       

%% Updates lineout plot
        function update_plots(hObject,eventdata,handles)
            try
            n=handles.list_idx;
            axes(handles.lineout_axes); cla(handles.lineout_axes);
            name = strsplit(handles.fnames{n},'_Ch'); name = name{1};
            if handles.calc_bool ==1
                switch handles.plot_idx
                    case 1 %velocity
                        plot(handles.lineout_time{n},handles.velocity{n},'.b');
                        xlabel('time (ns)');ylabel('velocity (km/s)');
                    case 2 %Pressure
                        plot(handles.lineout_time{n},handles.pressure{n},'.g');
                        xlabel('time (ns)');ylabel('pressure (GPa)'); 
                    case 3 %Flux
                        plot(handles.lineout_time{n},handles.flux{n},'.c');
                        xlabel('time (ns)');ylabel('flux (TJ/m-s'); 
                    case 4 %Fluence
                        plot(handles.lineout_time{n},handles.fluence{n},'.k');
                        xlabel('time (ns)');ylabel('fluence (kJ/m^2)');
                    case 5 %Displacement
                        plot(handles.lineout_time{n},handles.displacement{n},'.r');
                        xlabel('time (ns)');ylabel('displacement (um)');
                end
                update_timing_box(hObject,eventdata,handles);
            else
                cla(handles.lineout_axes);
                handles.current_message = 'No velocity calculated yet';
            end
            title(name);
            plot_int(hObject,eventdata,handles);
            linkaxes([handles.lineout_axes,handles.interferogram_axes],'x');
            error('Allfunctions:NoError','No Error');
            catch ME
                exception_handler(ME,hObject,eventdata,handles);
            end
%% Updates lineout plot w/t=0 set as first movement
%Note, I will update this later with a better time0 fitting mechanism (this
%one looks for 5RMS above noise), but for now this works pretty well.
 function update_plots_IM(hObject,eventdata,handles)
            try
            n=handles.list_idx;
            axes(handles.lineout_axes); cla(handles.lineout_axes);
            name = strsplit(handles.fnames{n},'_Ch'); name = name{1};
            if handles.calc_bool ==1
                new_ltime = handles.lineout_time{n}-handles.time{n}(handles.t0{n});
                switch handles.plot_idx
                    case 1 %velocity
                        plot(new_ltime,handles.velocity{n},'.b'); 
                        xlabel('time (ns)');ylabel('velocity (km/s)');
                    case 2 %Pressure
                        plot(new_ltime,handles.pressure{n},'.g');
                        xlabel('time (ns)');ylabel('pressure (GPa)');
                    case 3 %Flux
                        plot(new_ltime,handles.flux{n},'.c'); 
                        xlabel('time (ns)');ylabel('flux (TJ/m-s'); 
                    case 4 %Fluence
                        plot(new_ltime,handles.fluence{n},'.k');
                        xlabel('time (ns)');ylabel('fluence (kJ/m^2)'); 
                    case 5 %Displacement
                        plot(new_ltime,handles.displacement{n},'.r'); 
                        xlabel('time (ns)');ylabel('displacement (um)');
                end
                update_timing_box(hObject,eventdata,handles);
            else
                cla(handles.lineout_axes);
                handles.current_message = 'No velocity calculated yet';
            end
            title(name);
            plot_int_IM(hObject,eventdata,handles);
            linkaxes([handles.lineout_axes,handles.interferogram_axes],'x');
            error('Allfunctions:NoError','No Error');
            catch ME
                exception_handler(ME,hObject,eventdata,handles);
            end
            %%
 function int = plot_int_IM(hObject,eventdata,handles)
        n=handles.list_idx;
        amp = handles.amplitude{n};
        t0 = handles.time{n}(handles.t0{n});
        new_time = handles.time{n}-t0;
        axes(handles.interferogram_axes); hold off;
        for i = 1:4
            if i==3
                continue
            end
            plot(new_time,amp(:,i)); hold on;
            if (handles.method ==2 && handles.peaks_bool ==1)
                if i==4
                    xvals = handles.xyPeaks{n}{3}(:,1)-t0;yvals = handles.xyPeaks{n}{3}(:,2);
                else
                    xvals = handles.xyPeaks{n}{i}(:,1)-t0;yvals = handles.xyPeaks{n}{i}(:,2);
                end
                plot(xvals,yvals,'ok')
            end
        end
        xlabel('time (ns)');ylabel('volts (V)');
        if handles.peaks_bool ==0
            legend('Channel 1','Channel 2','Channel 3');
        else
            legend('Channel 1','','Channel 2','','Channel 3','');
        end
%% Deletes old data if new data is found
    function Delete_Old_Data(hObject,eventdata,handles)
        clear handles.velocity handles.lineout_time handles.xPeaks handles.yPeaks handles.xyPeaks
        handles.velocity = {};handles.lineout_time = {}; handles.xPeaks = {}; handles.yPeaks = {}; handles.xyPeaks = {};
        guidata(hObject,handles);
        
%% Update Timing Information Box

        function update_timing_box(hObject,eventdata,handles)
            n = handles.list_idx;
            %[~,min_idx] = min(handles.phase{n});
            min_idx = handles.min_idx{n};
            min_time = handles.time{n}(min_idx);
            output_text = sprintf('Oscilloscope offset = %s \n First Rise Time = %s \n Phase Minimum = %s',string(handles.scope_offset),...
                string(handles.time{n}(handles.t0{n})),string(min_time));
            set(handles.timing_info,'String',output_text);
%% aggregate error handler
    function error_out  = exception_handler(error_in,hObject,eventdata,handles)
        try
            line_num = error_in.stack.line;
        catch
            line_num = 0;
        end
        switch error_in.identifier
            case 'LoadData:NoData'
                output_text = 'No file loaded';
            case 'MATLAB:nonExistentField'
                if line_num ==144
                output_text = 'No Data to Run';
                else
                    output_text = sprintf('Unknown error: %s on line %d',...
                        string(error_in.message),line_num);
                end
            case 'Plotting:NoVelFile'
                output_text = 'Velocity not calculated';
                cla(handles.lineout_axes); 
            case 'LISTBOX:NOFILE'
                output_text = 'No file loaded';
            case 'MethodUpdate:Badmethod'
                output_text = 'Bad method number';
            case 'RunData:NoData'
                output_text = 'There is no data';
            case 'SaveFunction:NoSaveData'
                output_text = 'No Data to save';
            case 'SaveFunc:CancelInput'
                output_text = 'Save Cancelled';
            case 'Allfunctions:NoError'
                output_text = handles.current_message;
            otherwise
                output_text = sprintf('Unknown error: %s on line %d',string(error_in.message),line_num);
        end
        set(handles.text_prompt,'String',output_text)
 
        
        
        
   %% References and acknowledgements    
 %{
   Acknowledgements
   William Shaw and Alex Curtis for developing the original STFT analysis
   code (still mostly intact in this iteration);
   All work was performed under the advisement of Dana Dlott.
   
   References:
   [1] D. H. Dolan, S. C. Jones. THRIVE: a data reduction program for
        three-phase PDV/PDI and VISAR measurements
   [2] B. J. Jensen, D. B. Holtkamp, P. A. Rigg, D. H. Dolan. Accuracy
        limits and window corrections for photon Doppler velocimetry. Journal of
        Applied Physics, 101, 013523, 2007.
   [3] D. J. Chapman, D. E. Eakins, D. M. Williamson and W. G. Proud. Index of refraction
        measurements and window corrections for PMMA under shock compression.
        In M. L. Elert, W. T. Buttler, J. P. Borg, J. L. Jordan and T. J. Vogler, editors,
        Shock Compression Condens. Matter, volume 442, pages 442–445. AIP Conf. Proc.,
        Chicago, 2012. doi:10.1063/1.3686313
   [4] J. W. Forbes. Shock Wave Compression of Condensed Matter
   %}




