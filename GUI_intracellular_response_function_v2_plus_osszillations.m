function varargout = GUI_intracellular_response_function_v2_plus_osszillations(varargin)
% GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS MATLAB code for GUI_intracellular_response_function_v2_plus_osszillations.fig
%      GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS, by itself, creates a new GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS or raises the existing
%      singleton*.
%
%      H = GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS returns the handle to a new GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS or the handle to
%      the existing singleton*.
%
%      GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS.M with the given input arguments.
%
%      GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS('Property','Value',...) creates a new GUI_INTRACELLULAR_RESPONSE_FUNCTION_V2_PLUS_OSSZILLATIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_intracellular_response_function_v2_plus_osszillations_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_intracellular_response_function_v2_plus_osszillations_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the axes_response to help GUI_intracellular_response_function_v2_plus_osszillations

% Last Modified by GUIDE v2.5 29-Nov-2013 16:04:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_intracellular_response_function_v2_plus_osszillations_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_intracellular_response_function_v2_plus_osszillations_OutputFcn, ...
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


% --- Executes just before GUI_intracellular_response_function_v2_plus_osszillations is made visible.
function GUI_intracellular_response_function_v2_plus_osszillations_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_intracellular_response_function_v2_plus_osszillations (see VARARGIN)

% Choose default command line output for GUI_intracellular_response_function_v2_plus_osszillations
handles.output = hObject;

%set the default path
handles.dir='E:\Science\data\response_function\oocytes\2013-09-20\';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_intracellular_response_function_v2_plus_osszillations wait for user axes_response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_intracellular_response_function_v2_plus_osszillations_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_particle_diameter_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to edit_particle_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_particle_diameter as text
%        str2double(get(hObject,'String')) returns contents of edit_particle_diameter as a double
handles=display_current_data(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_particle_diameter_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to edit_particle_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_response_prefactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_response_prefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_response_prefactor as text
%        str2double(get(hObject,'String')) returns contents of edit_response_prefactor as a double
handles=display_current_data(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_response_prefactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_response_prefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_folder_Callback(hObject, eventdata, handles)
% hObject    handle to load_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.dir=uigetdir(handles.dir,'Give me the base directory of the response data');
folders=dir([handles.dir]);
folders(1:2)=[];
folders(find([folders.isdir]==0))=[];
handles.folders=folders;
handles.current_folder=1;
set(handles.edit_dataset,'String','1');
handles.current_path=folders(handles.current_folder);
set(handles.dataset_slider,'Value',1);
set(handles.dataset_slider,'Max',length(folders));
set(handles.dataset_slider,'Min',1);
set(handles.dataset_slider,'SliderStep',[1,1]/(length(folders)-1));

%Here I make sure that the collected response strcuture is empty
handles.corrected_response=[];


handles=load_current_folder(handles);

guidata(hObject, handles);

%--------------------
function handles=load_current_folder(handles)
%This function loads the data in the current folder

handles.current_folder=str2num(get(handles.edit_dataset,'String'));
folders=handles.folders;
%now I read all the response folders and populate the listbox
folders_p=dir([handles.dir,filesep,folders(handles.current_folder).name,filesep,'multiple_run_series*']);
set(handles.listbox_passive,'String',{folders_p.name});

handles.folders_p=folders_p;
%now I read all the fluctuations folders and populate the listbox
folders_a=dir([handles.dir,filesep,folders(handles.current_folder).name,filesep,'response_function*']);
set(handles.listbox_active,'String',{folders_a.name});

handles.folders_a=folders_a;
%Finally I define the first datasets and load these as current data
set(handles.listbox_passive,'Value',1);
set(handles.listbox_active,'Value',1);
handles=load_current_dataset(handles);


function handles=load_current_dataset(handles)
%Here I will take the datasets marked in the listboxes and load the actual
%data. After that I will directly display them

%load the active response
folders=handles.folders;
folders_a=handles.folders_a;
folders_p=handles.folders_p;
act_data_a=get(handles.listbox_active,'Value');
act_data_p=get(handles.listbox_passive,'Value');

load_active=[handles.dir,filesep,folders(handles.current_folder).name,filesep,folders_a(act_data_a).name];
try
    load_passive=[handles.dir,filesep,folders(handles.current_folder).name,filesep,folders_p(act_data_p).name];
end
[f,alphax,alphay,trap_stiff,slopes,act_trap]=process_response_function_folder(load_active);

%here I will get the characteristics of the current active data
handles.active_files=dir([load_active,filesep,'*deformation_response.mat']);
%scale the screen frequency slider
set(handles.slider_freq,'Value',1);
set(handles.slider_freq,'Max',length(handles.active_files));
set(handles.slider_freq,'Min',1);
set(handles.slider_freq,'SliderStep',[1,1]/(length(handles.active_files)-1));


handles.response_f=f;
handles.response_x=alphax;
handles.response_y=alphay;
handles.trap_stiff=trap_stiff;
handles.slopes=slopes;
handles.act_trap=act_trap;


%now I will either load the fluctuations that was saved with the active
%microrheology, or I will take the fluctuations from the passive
%microrheology folder
if isempty(handles.folders_p)
    load([load_active,filesep,'histogram_results.mat']);
    handles.psd_f=fl(1,2:end);
    handles.psd_x=pxl;
    handles.psd_y=pyl;
else
    [fluct,fl,psdxy,psdxy_std]=process_fluctuation_folder_v3(load_passive); %#ok<NASGU,ASGLU>
    handles.psd_f=fl(1,:);
    handles.psd_x=psdxy(1,:);
    handles.psd_y=psdxy(2,:);
end


handles.particle_diameter=str2num(get(handles.edit_particle_diameter,'String'))*1e-6;
handles.response_prefactor=str2num(get(handles.edit_response_prefactor,'String'));

handles=display_current_data(handles);


%--------------------
function handles=display_current_data(handles)
%This function displays the data in the current folder

handles.particle_diameter=str2num(get(handles.edit_particle_diameter,'String'))*1e-6;
handles.response_prefactor=str2num(get(handles.edit_response_prefactor,'String'));

if not(isfield(handles , 'response_function_dirrection'))
    handles.response_function_dirrection = 'X';
end
    
if handles.response_function_dirrection == 'X'
    response_xory =handles.response_x;
elseif handles.response_function_dirrection == 'Y'
    response_xory =handles.response_y;
else
    throw('unknown direction')
end
psd_x=handles.psd_x(handles.act_trap,:);

axes(handles.axes_response)


eta=1./(3*pi*handles.particle_diameter*2*pi*handles.response_f.*handles.response_prefactor .* abs(imag(response_xory)));
loglog(handles.response_f,eta);
hold on
G=1./(6*pi*handles.particle_diameter*handles.response_prefactor .* ((response_xory)));
loglog(handles.response_f,abs(real(G)),'+r-');
loglog(handles.response_f,abs(imag(G)),'or-');
loglog([100:1000],[100:1000].^.75*.1);

legend('viscosity','real(G)','img(G)','x^3/4')

xlabel('f in [Hz]')
ylabel('G in Pa (o dissipative, + elastic)')




freq=length(G)-round(get(handles.slider_freq,'Value'))+1;
%now I load the frequency selected
%load the active response
folders=handles.folders;
folders_a=handles.folders_a;
act_data_a=get(handles.listbox_active,'Value');
load_active=[handles.dir,filesep,folders(handles.current_folder).name,filesep,folders_a(act_data_a).name];
load([load_active,filesep,handles.active_files(freq).name]);
data=squeeze(data);

%I quickly plot a line a the current frequency
plot([f f],[0.001 100],'k')

hold off

axes(handles.oss);
%first I plot the normaized force and position
xy_k=-xy_k;
x=-data(1,:)/cal(5)*1e-6+1./(xy_slope(1).*1e6).*data(4,:)./data(8,:);
Fx=xy_k(1)./(xy_slope(1)*1e6).*data(4,:)./data(8,:);

%lets fit a sine to the force
temp=[0:length(x)-1]/s_eff;
f_sin_f=fittype('a*sin(2*pi*f*t+phase)','independent','t','problem','f')
[ff,goff]=fit(temp',Fx'-mean(Fx),f_sin_f,'problem',f)
phase=ff.phase


f_sin_p=fittype('a*sin(2*pi*f*t+phase)+b*cos(2*pi*f*t+phase)','independent','t','problem',{'f','phase'});
[fp,gofp]=fit(temp',x'-mean(x),f_sin_p,'problem',{f,phase})
clear i;
response_recal=fp.a/ff.a+i*(fp.b/ff.a)

%I now calculate the length of the first 8 periods
l_l=8*s_eff/f;
%plotyy([1:length(data)]/s_eff,x,[1:length(data)]/s_eff,Fx);
plotyy([1:l_l],x(1:l_l),[1:l_l],Fx(1:l_l));



%now I get the fourier transform of this data ad plot the force and the
%position
p=length(x);
FX=fft(x);
FF=fft(Fx);
FX=FX(1:p/2+1);
FF=FF(1:p/2+1);
alpha_x=fft(x)./fft(Fx);
alpha_x=alpha_x(1:p/2+1);
fr=s_eff/p*([0:p/2]);


axes(handles.fft);
loglog(fr,real(FX),'b-+')
hold on
loglog(fr,imag(FX),'b-o')
loglog(fr,real(FF),'g-+')
loglog(fr,imag(FF),'g-o')
xlim([min(fr) max(fr)]);
ylim([1e-12 1e-2])
hold off

%now I get teh absolute value and the phase of the actualt frequency
[a,b]=min(abs(f-fr));
cx=FX(b);
cf=FF(b);
fx_p=atan2(imag(cf),real(cf))
fx_a = abs(cf)
x_p=atan2(imag(cx),real(cx))
x_a = abs(cx)





save_data.particle_diameter=handles.particle_diameter;
save_data.response_prefactor=handles.response_prefactor;
save_data.response_f=handles.response_f;
save_data.corr_active_response=abs(handles.response_prefactor * imag(response_xory));
save_data.psd_f=handles.psd_f;
save_data.passive_response=abs(pi.*handles.psd_f./4e-21.*psd_x);
save_data.G=G;
save_data.eta=eta;
save_data.base_path=handles.dir;
save_data.folder_path=folders(handles.current_folder).name;

%and now store this structure in the handles.
handles.save_data=save_data;







function edit_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dataset as text
%        str2double(get(hObject,'String')) returns contents of edit_dataset as a double
set(handles.dataset_slider,'Value',str2num(get(handles.edit_dataset,'String')));
handles=load_current_folder(handles);
handles=display_current_data(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_current.
function save_current_Callback(hObject, eventdata, handles)
% hObject    handle to save_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Here I will save the currently discplayed data in a file that is in the
%current path.
save_data=handles.save_data;
try 
    corrected_response=handles.corrected_response;
    i=length(corrected_response);
    corrected_response(i+1)=save_data;
catch
    i=0;
    corrected_response=save_data;
end

handles.corrected_response=corrected_response;
guidata(hObject, handles);




% --- Executes on slider movement.
function dataset_slider_Callback(hObject, eventdata, handles)
% hObject    handle to dataset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v=get(handles.dataset_slider,'Value');
set(handles.edit_dataset,'String',num2str(v));
handles=load_current_folder(handles);
handles=display_current_data(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dataset_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in listbox_passive.
function listbox_passive_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_passive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_passive contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_passive
handles=load_current_dataset(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox_passive_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_passive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_active.
function listbox_active_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_active contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_active
handles=load_current_dataset(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox_active_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in save_all.
function save_all_Callback(hObject, eventdata, handles)
% hObject    handle to save_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

corrected_response=handles.corrected_response;
folders=handles.folders;
[s_file,s_path]=uiputfile([handles.dir,filesep,'*.mat'],'Where should I store the data');
save([s_path,filesep,s_file],'corrected_response');


% --- Executes when selected object is changed in XorY.
function XorY_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in XorY 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
    switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
        case 'radiobuttonX'
        handles.response_function_dirrection = 'X';
        case 'radiobuttonY'
        handles.response_function_dirrection = 'Y';
        % Code for when radiobutton2 is selected.
        otherwise
        throw('got unknown radio button !') 
    end
    guidata(hObject, handles);


% --- Executes on slider movement.
function slider_freq_Callback(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles=display_current_data(handles)

 guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
