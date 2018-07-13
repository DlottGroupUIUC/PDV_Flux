function varargout = PDV_flux_3_TIM(varargin)
% PDV_FLUX_3_TIM MATLAB code for PDV_flux_3_TIM.fig
%      PDV_FLUX_3_TIM, by itself, creates a new PDV_FLUX_3_TIM or raises the existing
%      singleton*.
%
%      H = PDV_FLUX_3_TIM returns the handle to a new PDV_FLUX_3_TIM or the handle to
%      the existing singleton*.
%
%      PDV_FLUX_3_TIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PDV_FLUX_3_TIM.M with the given input arguments.
%
%      PDV_FLUX_3_TIM('Property','Value',...) creates a new PDV_FLUX_3_TIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PDV_flux_3_TIM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PDV_flux_3_TIM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PDV_flux_3_TIM

% Last Modified by GUIDE v2.5 09-Jul-2018 18:05:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PDV_flux_3_TIM_OpeningFcn, ...
                   'gui_OutputFcn',  @PDV_flux_3_TIM_OutputFcn, ...
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


% --- Executes just before PDV_flux_3_TIM is made visible.
function PDV_flux_3_TIM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PDV_flux_3_TIM (see VARARGIN)

% Choose default command line output for PDV_flux_3_TIM
handles.output = hObject;
h = findobj('Name','PDV_flux_3');
if ~isempty(h)
    main_GUI_data = guidata(h);
end
handles.main_GUI_data = main_GUI_data;
gd = handles.main_GUI_data; n = gd.list_idx;
handles.tr = [gd.time{n}(gd.t0{n})-20,gd.time{n}(gd.t0{n})+20];
handles.offset = gd.time_offset;
% Update handles structure
guidata(hObject, handles);
plot_int_tim(hObject,eventdata,handles);
% UIWAIT makes PDV_flux_3_TIM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PDV_flux_3_TIM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_button.
function start_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x1,~] = ginput(1);
%[x2,~] = ginput(1);
h = handles.main_GUI_data; n = h.list_idx;
amp = h.camp{n};
[~,t_center] = min(abs(h.time{n}-x1));
handles.t_center = h.time{n}(t_center);
ti = t_center-250;tf = t_center+50;
handles.tr = [h.time{n}(ti),h.time{n}(tf)];
guidata(hObject,handles);
plot_int_tim(hObject,eventdata,handles);

[handles.t,handles.intensity] = STFT(ti,tf,h.time{n},amp,handles.offset{n});
[logitC,dev] = glmfit(handles.t,handles.intensity,'binomial','logit');
handles.logitFit = glmval(logitC,handles.t,'logit');
plot_intensity(hObject,eventdata,handles);


% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(PDV_flux_3_TIM);
%%Begin Subroutine Section

%% Find TIM based on hamming window spectrogram [5]
    function [t,intensity_fit] = STFT(ti,tf,time,amp,offset)
        ntime = time.*1e-9;
        sample_frequency = 1./(ntime(2)-ntime(1));
        r = round(sample_frequency.*(2*10^-9));
        toff = time(ti);ampi = amp(ti:tf);
        k=1;
        for i=1:4
            if i==3
                continue
            end
        [STFT(:,:,k),f,t] = spectrogram(ampi,hamming(r),r-1,64*r,sample_frequency);
        k=  k+1;
        end
        STFT = (abs(STFT(:,:,1))+abs(STFT(:,:,2))+abs(STFT(:,:,3)))./3;
        t = t.*1e9+toff;
        [mx locs]=max(STFT,[],1);
        intensity=f(locs);

        % Fit the FFT at each time step to better resolve the velocity. I use a
        % polynomial since this is much, much less computationally expensive then a
        % gaussian fit. 
        intensity_fit = intensity;
        %{
        for i=1:length(intensity)
            if (locs(i)+2)<length(f)&&(locs(i)-2)>1
                p = polyfit(f((locs(i)-2):(locs(i)+2)),STFT((locs(i)-2):(locs(i)+2),i),2);
                peakPosition = -p(2)./(p(1)*2);
                intensity_fit(i) = peakPosition;
            else
                intensity_fit(i) = 0;
            end
        end
        %}
        intensity_fit = intensity_fit./max(intensity_fit);

        

%% plots interferogram
    function int = plot_int_tim(hObject,eventdata,handles)
        h=handles.main_GUI_data;  n=h.list_idx;
        name = strsplit(h.fnames{n},'_Ch');
        name = name{1};
        amp = h.amplitude{n};
        axes(handles.main_axes); hold off;
        
        plot(h.time{n},amp(:,2));
        xlabel('time (ns)');ylabel('volts (V)'); title(sprintf('%s: %s',name,'TIM'));
        xlim(handles.tr);
        legend('Channel 2');
        
    function int = plot_intensity(hObject,eventdata,handles)
        h=handles.main_GUI_data;  n=h.list_idx;
        name = strsplit(h.fnames{n},'_Ch');
        name = name{1};
        axes(handles.main_axes); hold off;
        
        plot(handles.t,handles.intensity,'k'); hold on;
        plot(handles.t,handles.logitFit,'r');
        xlabel('time (ns)');ylabel('intensity (normalized)'); title(sprintf('%s: %s',name,'TIM'));
        tr_plot = [handles.t_center-1,handles.t_center+1];
        xlim(tr_plot);
        legend('STFT 2ns window','Logit Fit');
        
     function int = plot_int_tim_2(hObject,eventdata,handles)
        h=handles.main_GUI_data;  n=h.list_idx;
        name = strsplit(h.fnames{n},'_Ch');
        name = name{1};
        amp = h.amplitude{n};
        axes(handles.main_axes); hold off;
        
        plot(h.time{n},amp(:,2)); hold on;
        plot(h.time{n}(handles.ti:handles.tf),handles.logitFit,'r');
        xlabel('time (ns)');ylabel('volts (V)'); title(sprintf('%s: %s',name,'TIM'));
        tr_plot = [handles.t_center-3,handles.t_center+3];
        xlim(tr_plot);
        legend('STFT 2ns window','Logit Fit');
