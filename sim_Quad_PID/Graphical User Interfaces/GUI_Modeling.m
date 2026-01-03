function varargout = GUI_Modeling(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Modeling_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Modeling_OutputFcn, ...
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


% --- Executes just before GUI_Modeling is made visible.
function GUI_Modeling_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% --- کد نمایش لوگو ---
try
    % نام فایل لوگوی خود را اینجا وارد کنید
    logo = imread('your_logo.JPG'); 
    % نمایش در محوری که تگ آن logoAxes است
    axes(handles.logoAxes);
    imshow(logo);
    axis(handles.logoAxes, 'off');
catch
    disp('فایل لوگو پیدا نشد. لطفاً نام و مسیر فایل را چک کنید.');
end

% --- کد اصلاح شده نمایش تصویر وسط ---
% با استفاده از تگ صحیحی که پیدا کردید: 'graphics'
imshow('Motors_Graphic.jpg', 'Parent', handles.graphics); 

set(handles.motor_m_unit, 'String', 'g');
set(handles.ESC_m_unit, 'String', 'g');
set(handles.HUB_m_unit, 'String', 'g');
set(handles.Arms_m_unit, 'String', 'g');
set(handles.motor_dm_unit, 'String', 'cm');
set(handles.motor_h_unit, 'String', 'cm');
set(handles.motor_r_unit, 'String', 'cm');
set(handles.ESC_a_unit, 'String', 'cm');
set(handles.ESC_b_unit, 'String', 'cm');
set(handles.ESC_ds_unit, 'String', 'cm');
set(handles.HUB_r_unit, 'String', 'cm');
set(handles.HUB_H_unit, 'String', 'cm');
set(handles.Arms_r_unit, 'String', 'cm');
set(handles.Arms_L_unit, 'String', 'cm');
set(handles.Arms_da_unit, 'String', 'cm');
movegui('center');
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Modeling_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function massTotal_edit_Callback(hObject, eventdata, handles)
function massTotal_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ct_edit_Callback(hObject, eventdata, handles)
function ct_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cq_edit_Callback(hObject, eventdata, handles)
function cq_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function motor_m_edit_Callback(hObject, eventdata, handles)
function motor_m_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function motor_dm_edit_Callback(hObject, eventdata, handles)
function motor_dm_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function motor_h_edit_Callback(hObject, eventdata, handles)
function motor_h_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function motor_r_edit_Callback(hObject, eventdata, handles)
function motor_r_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ESC_m_edit_Callback(hObject, eventdata, handles)
function ESC_m_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ESC_a_edit_Callback(hObject, eventdata, handles)
function ESC_a_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ESC_b_edit_Callback(hObject, eventdata, handles)
function ESC_b_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ESC_ds_edit_Callback(hObject, eventdata, handles)
function ESC_ds_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function arms_m_edit_Callback(hObject, eventdata, handles)
function arms_m_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function arms_r_edit_Callback(hObject, eventdata, handles)
function arms_r_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function arms_L_edit_Callback(hObject, eventdata, handles)
function arms_L_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HUB_m_edit_Callback(hObject, eventdata, handles)
function HUB_m_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HUB_r_edit_Callback(hObject, eventdata, handles)
function HUB_r_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HUB_H_edit_Callback(hObject, eventdata, handles)
function HUB_H_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes when selected object is changed in pic_display_panel.
function pic_display_panel_SelectionChangeFcn(hObject, eventdata, handles)
% --- کد اصلاح شده برای دکمه‌های رادیویی ---
% با استفاده از تگ صحیحی که پیدا کردید: 'graphics'
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object
    case 'motors_button'
        imshow('Motors_Graphic.jpg', 'Parent', handles.graphics);
    case 'ESC_button'
        imshow('ESC_Graphic.jpg', 'Parent', handles.graphics);
    case 'HUB_button'
        imshow('HUB_Graphic.jpg', 'Parent', handles.graphics);
    case 'arms_button'
        imshow('ARMS_Graphic.jpg', 'Parent', handles.graphics);
end
function jx_Callback(hObject, eventdata, handles)
function jx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function jy_Callback(hObject, eventdata, handles)
function jy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function jz_Callback(hObject, eventdata, handles)
function jz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in save_plus.
function save_plus_Callback(hObject, eventdata, handles)
massUnits = get(handles.motor_m_unit, 'String');
if massUnits == 'g',
    motor_m_g = str2num(get(handles.motor_m_edit, 'String'));
    motor_dm_cm = str2num(get(handles.motor_dm_edit, 'String'));
    motor_h_cm = str2num(get(handles.motor_h_edit, 'String'));
    motor_r_cm = str2num(get(handles.motor_r_edit, 'String'));
    ESC_m_g = str2num(get(handles.ESC_m_edit, 'String'));
    ESC_a_cm = str2num(get(handles.ESC_a_edit, 'String'));
    ESC_b_cm = str2num(get(handles.ESC_b_edit, 'String'));
    ESC_ds_cm = str2num(get(handles.ESC_ds_edit, 'String'));
    HUB_m_g = str2num(get(handles.HUB_m_edit, 'String'));
    HUB_r_cm = str2num(get(handles.HUB_r_edit, 'String'));
    HUB_H_cm = str2num(get(handles.HUB_H_edit, 'String'));
    arms_m_g = str2num(get(handles.arms_m_edit, 'String'));
    arms_r_cm = str2num(get(handles.arms_r_edit, 'String'));
    arms_L_cm = str2num(get(handles.arms_L_edit, 'String'));
    arms_da_cm = str2num(get(handles.arms_da_edit, 'String'));
    motor_m = motor_m_g/1000;
    motor_dm = motor_dm_cm/100;
    motor_h = motor_h_cm/100;
    motor_r = motor_r_cm/100;
    ESC_m = ESC_m_g/1000;
    ESC_a = ESC_a_cm/100;
    ESC_b = ESC_b_cm/100;
    ESC_ds = ESC_ds_cm/100;
    HUB_m = HUB_m_g/1000;
    HUB_r = HUB_r_cm/100;
    HUB_H = HUB_H_cm/100;
    arms_m = arms_m_g/1000;
    arms_r = arms_r_cm/100;
    arms_L = arms_L_cm/100;
    arms_da = arms_da_cm/100;
    d = motor_dm;
else 
    motor_m_oz = str2num(get(handles.motor_m_edit, 'String'));
    motor_dm_in = str2num(get(handles.motor_dm_edit, 'String'));
    motor_h_in = str2num(get(handles.motor_h_edit, 'String'));
    motor_r_in = str2num(get(handles.motor_r_edit, 'String'));
    ESC_m_oz = str2num(get(handles.ESC_m_edit, 'String'));
    ESC_a_in = str2num(get(handles.ESC_a_edit, 'String'));
    ESC_b_in = str2num(get(handles.ESC_b_edit, 'String'));
    ESC_ds_in = str2num(get(handles.ESC_ds_edit, 'String'));
    HUB_m_oz = str2num(get(handles.HUB_m_edit, 'String'));
    HUB_r_in = str2num(get(handles.HUB_r_edit, 'String'));
    HUB_H_in = str2num(get(handles.HUB_H_edit, 'String'));
    arms_m_oz = str2num(get(handles.arms_m_edit, 'String'));
    arms_r_in = str2num(get(handles.arms_r_edit, 'String'));
    arms_L_in = str2num(get(handles.arms_L_edit, 'String'));
    arms_da_in = str2num(get(handles.arms_da_edit, 'String'));
    motor_m = motor_m_oz/35.274;
    motor_dm = motor_dm_in/39.3701;
    motor_h = motor_h_in/39.3701;
    motor_r = motor_r_in/39.3701;
    ESC_m = ESC_m_oz/35.274;
    ESC_a = ESC_a_in/39.3701;
    ESC_b = ESC_b_in/39.3701;
    ESC_ds = ESC_ds_in/39.3701;
    HUB_m = HUB_m_oz/35.274;
    HUB_r = HUB_r_in/39.3701;
    HUB_H = HUB_H_in/39.3701;
    arms_m = arms_m_oz/35.274;
    arms_r = arms_r_in/39.3701;
    arms_L = arms_L_in/39.3701;
    arms_da = arms_da_in/39.3701;
    d = motor_dm;
end
g = 9.81;
mass = str2num(get(handles.massTotal_edit, 'String'));
T = str2num(get(handles.tConstant, 'String'));
minThr = str2num(get(handles.minThr, 'String'));
cr = str2num(get(handles.cr_edit, 'String'));
b = str2num(get(handles.b_edit, 'String'));
ct = str2num(get(handles.ct_edit, 'String'));
cq = str2num(get(handles.cq_edit, 'String'));
Jx = str2num(get(handles.jx, 'String'));
Jy = str2num(get(handles.jy, 'String'));
Jz = str2num(get(handles.jz, 'String'));
Jb = [Jx 0 0; 0 Jy 0; 0 0 Jz];
Jbinv = [1/Jx 0 0; 0 1/Jy 0; 0 0 1/Jz];
dctcq = [0 d*ct 0 -d*ct; -d*ct 0 d*ct 0; -cq cq -cq cq];
plusConfig = 1;
mRC = (motor_m)*(0.527);
Jm = ((mRC)*(motor_r)^2)/2;
quadModel = struct('g',(g),'d',(d),'mass',(mass),'ct',(ct),'cq',(cq),...
    'Jx',(Jx),'Jy',(Jy),'Jz',(Jz),'Jm',(Jm),'Jb',(Jb),'Jbinv',(Jbinv),'dctcq',(dctcq),...
    'motor_m',(motor_m),'motor_dm',(motor_dm),'motor_h',(motor_h),'motor_r',(motor_r),...
    'ESC_m',(ESC_m),'ESC_a',(ESC_a),'ESC_b',(ESC_b),'ESC_ds',(ESC_ds),...
    'HUB_m',(HUB_m),'HUB_r',(HUB_r),'HUB_H',(HUB_H),...
    'arms_m',(arms_m),'arms_r',(arms_r),'arms_L',(arms_L), 'arms_da',(arms_da),'T',(T),'minThr',(minThr),...
    'cr',(cr),'b',(b), 'plusConfig',(plusConfig));
uisave('quadModel','quadModel_+');
guidata(hObject, handles);
% --- Executes on button press in save_X.
function save_X_Callback(hObject, eventdata, handles)
massUnits = get(handles.motor_m_unit, 'String');
if massUnits == 'g',
    motor_m_g = str2num(get(handles.motor_m_edit, 'String'));
    motor_dm_cm = str2num(get(handles.motor_dm_edit, 'String'));
    motor_h_cm = str2num(get(handles.motor_h_edit, 'String'));
    motor_r_cm = str2num(get(handles.motor_r_edit, 'String'));
    ESC_m_g = str2num(get(handles.ESC_m_edit, 'String'));
    ESC_a_cm = str2num(get(handles.ESC_a_edit, 'String'));
    ESC_b_cm = str2num(get(handles.ESC_b_edit, 'String'));
    ESC_ds_cm = str2num(get(handles.ESC_ds_edit, 'String'));
    HUB_m_g = str2num(get(handles.HUB_m_edit, 'String'));
    HUB_r_cm = str2num(get(handles.HUB_r_edit, 'String'));
    HUB_H_cm = str2num(get(handles.HUB_H_edit, 'String'));
    arms_m_g = str2num(get(handles.arms_m_edit, 'String'));
    arms_r_cm = str2num(get(handles.arms_r_edit, 'String'));
    arms_L_cm = str2num(get(handles.arms_L_edit, 'String'));
    arms_da_cm = str2num(get(handles.arms_da_edit, 'String'));
    motor_m = motor_m_g/1000;
    motor_dm = motor_dm_cm/100;
    motor_h = motor_h_cm/100;
    motor_r = motor_r_cm/100;
    ESC_m = ESC_m_g/1000;
    ESC_a = ESC_a_cm/100;
    ESC_b = ESC_b_cm/100;
    ESC_ds = ESC_ds_cm/100;
    HUB_m = HUB_m_g/1000;
    HUB_r = HUB_r_cm/100;
    HUB_H = HUB_H_cm/100;
    arms_m = arms_m_g/1000;
    arms_r = arms_r_cm/100;
    arms_L = arms_L_cm/100;
    arms_da = arms_da_cm/100;
    d = motor_dm;
else 
    motor_m_oz = str2num(get(handles.motor_m_edit, 'String'));
    motor_dm_in = str2num(get(handles.motor_dm_edit, 'String'));
    motor_h_in = str2num(get(handles.motor_h_edit, 'String'));
    motor_r_in = str2num(get(handles.motor_r_edit, 'String'));
    ESC_m_oz = str2num(get(handles.ESC_m_edit, 'String'));
    ESC_a_in = str2num(get(handles.ESC_a_edit, 'String'));
    ESC_b_in = str2num(get(handles.ESC_b_edit, 'String'));
    ESC_ds_in = str2num(get(handles.ESC_ds_edit, 'String'));
    HUB_m_oz = str2num(get(handles.HUB_m_edit, 'String'));
    HUB_r_in = str2num(get(handles.HUB_r_edit, 'String'));
    HUB_H_in = str2num(get(handles.HUB_H_edit, 'String'));
    arms_m_oz = str2num(get(handles.arms_m_edit, 'String'));
    arms_r_in = str2num(get(handles.arms_r_edit, 'String'));
    arms_L_in = str2num(get(handles.arms_L_edit, 'String'));
    arms_da_in = str2num(get(handles.arms_da_edit, 'String'));
    motor_m = motor_m_oz/35.274;
    motor_dm = motor_dm_in/39.3701;
    motor_h = motor_h_in/39.3701;
    motor_r = motor_r_in/39.3701;
    ESC_m = ESC_m_oz/35.274;
    ESC_a = ESC_a_in/39.3701;
    ESC_b = ESC_b_in/39.3701;
    ESC_ds = ESC_ds_in/39.3701;
    HUB_m = HUB_m_oz/35.274;
    HUB_r = HUB_r_in/39.3701;
    HUB_H = HUB_H_in/39.3701;
    arms_m = arms_m_oz/35.274;
    arms_r = arms_r_in/39.3701;
    arms_L = arms_L_in/39.3701;
    arms_da = arms_da_in/39.3701;
    d = motor_dm;
end
g = 9.81;
mass = str2num(get(handles.massTotal_edit, 'String'));
T = str2num(get(handles.tConstant, 'String'));
minThr = str2num(get(handles.minThr, 'String'));
cr = str2num(get(handles.cr_edit, 'String'));
b = str2num(get(handles.b_edit, 'String'));
ct = str2num(get(handles.ct_edit, 'String'));
cq = str2num(get(handles.cq_edit, 'String'));
Jx = str2num(get(handles.jx, 'String'));
Jy = str2num(get(handles.jy, 'String'));
Jz = str2num(get(handles.jz, 'String'));
Jb = [Jx 0 0; 0 Jy 0; 0 0 Jz];
Jbinv = [1/Jx 0 0; 0 1/Jy 0; 0 0 1/Jz];
d45 = d*(sqrt(2)/2);
dctcq = [-d45*ct d45*ct d45*ct -d45*ct; -d45*ct -d45*ct d45*ct d45*ct; -cq cq -cq cq];
plusConfig = 0;
mRC = (motor_m)*(0.527);
Jm = ((mRC)*(motor_r)^2)/2;
quadModel = struct('g',(g),'d',(d),'mass',(mass),'ct',(ct),'cq',(cq),...
    'Jx',(Jx),'Jy',(Jy),'Jz',(Jz),'Jm',(Jm),'Jb',(Jb),'Jbinv',(Jbinv),'dctcq',(dctcq),...
    'motor_m',(motor_m),'motor_dm',(motor_dm),'motor_h',(motor_h),'motor_r',(motor_r),...
    'ESC_m',(ESC_m),'ESC_a',(ESC_a),'ESC_b',(ESC_b),'ESC_ds',(ESC_ds),...
    'HUB_m',(HUB_m),'HUB_r',(HUB_r),'HUB_H',(HUB_H),...
    'arms_m',(arms_m),'arms_r',(arms_r),'arms_L',(arms_L),'arms_da',(arms_da),'T',(T),'minThr',(minThr),...
    'cr',(cr),'b',(b),'plusConfig',(plusConfig));
uisave('quadModel','quadModel_X');
guidata(hObject, handles);
% --- Executes on button press in clear_button.
function clear_button_Callback(hObject, eventdata, handles)
set(handles.jx,'String','');
set(handles.jy,'String','');
set(handles.jz,'String','');
set(handles.massTotal_edit,'String','');
set(handles.ct_edit,'String','');
set(handles.cq_edit,'String','');
set(handles.motor_m_edit,'String','');
set(handles.motor_dm_edit,'String','');
set(handles.motor_h_edit,'String','');
set(handles.motor_r_edit,'String','');
set(handles.ESC_m_edit,'String','');
set(handles.ESC_a_edit,'String','');
set(handles.ESC_b_edit,'String','');
set(handles.ESC_ds_edit,'String','');
set(handles.HUB_m_edit,'String','');
set(handles.HUB_r_edit,'String','');
set(handles.HUB_H_edit,'String','');
set(handles.arms_m_edit,'String','');
set(handles.arms_r_edit,'String','');
set(handles.arms_L_edit,'String','');
set(handles.arms_da_edit,'String','');
set(handles.tConstant,'String','');
set(handles.minThr,'String','');
set(handles.cr_edit,'String','');
set(handles.b_edit,'String','');
% --- Executes on button press in calulateButton.
function calulateButton_Callback(hObject, eventdata, handles)
massUnits = get(handles.motor_m_unit, 'String');
if massUnits == 'g',
    mm_g = str2num(get(handles.motor_m_edit,'String'));
    dm_cm = str2num(get(handles.motor_dm_edit,'String'));
    hm_cm = str2num(get(handles.motor_h_edit,'String'));
    rm_cm = str2num(get(handles.motor_r_edit,'String')); 
    mm = mm_g/1000;
    dm = dm_cm/100;
    hm = hm_cm/100;
    rm = rm_cm/100;
    Jx1 = (mm) * ((rm)^2) + (4/3)*(mm) * ((hm)^2) + 2*(mm)*((dm)^2);
    Jy1 = (mm) * ((rm)^2) + (4/3)*(mm) * ((hm)^2) + 2*(mm)*((dm)^2);
    Jz1 = 2*(mm)*((rm)^2) + (4)*(mm)*((dm)^2);
    ms_g = str2num(get(handles.ESC_m_edit,'String'));
    as_cm = str2num(get(handles.ESC_a_edit,'String')); 
    bs_cm = str2num(get(handles.ESC_b_edit,'String')); 
    ds_cm = str2num(get(handles.ESC_ds_edit,'String')); 
    ms = ms_g/1000;
    as = as_cm/100;
    bs = bs_cm/100;
    ds = ds_cm/100;
    Jx2 = ((ms*(as^2))/6) + ((ms*(bs^2))/6) + (2*(ms)*(ds^2));
    Jy2 = ((ms*(as^2))/6) + ((ms*(bs^2))/6) + (2*(ms)*(ds^2));
    Jz2 = ((ms*(as^2 + bs^2))/3) + 4*(ms)*(ds^2);
    mh_g = str2num(get(handles.HUB_m_edit,'String')); 
    rh_cm = str2num(get(handles.HUB_r_edit,'String')); 
    Hh_cm = str2num(get(handles.HUB_H_edit,'String')); 
    mh = mh_g/1000;
    rh = rh_cm/100;
    Hh = Hh_cm/100;
    Jx3 = (1/4)*((mh))*((rh)^2) + (1/12)*((mh))*((Hh)^2);
    Jy3 = (1/4)*((mh))*((rh)^2) + (1/12)*((mh))*((Hh)^2);
    Jz3 = (1/2)*((mh))*((rh)^2);
    ma_g = str2num(get(handles.arms_m_edit,'String'));
    ra_cm = str2num(get(handles.arms_r_edit,'String'));
    La_cm = str2num(get(handles.arms_L_edit,'String'));
    da_cm = str2num(get(handles.arms_da_edit,'String'));
    ma = ma_g/1000;
    ra = ra_cm/100;
    La = La_cm/100;
    da = da_cm/100;
    Jx4 = ((3/2)*(ma)*(ra^2)) + ((2/3)*(ma)*(La^2)) + (2)*(ma)*(da^2);
    Jy4 = ((3/2)*(ma)*(ra^2)) + ((2/3)*(ma)*(La^2)) + (2)*(ma)*(da^2);
    Jz4 = ((ma)*((ra)^2)) + (4/3)*((ma) * ((La)^2)) + (4)*(ma)*(da^2);
    totalX = Jx1 + Jx2 + Jx3 + Jx4;
    totalY = Jy1 + Jy2 + Jy3 + Jy4;
    totalZ = Jz1 + Jz2 + Jz3 + Jz4;
    a = num2str(totalX);
    b = num2str(totalY);
    c = num2str(totalZ);
    set(handles.jx,'String',a);
    set(handles.jy,'String',b);
    set(handles.jz,'String',c);
    massTotal = 4*mm + 4*ms + mh + 4*ma;
    set(handles.massTotal_edit,'String',massTotal);
else 
    mm_oz = str2num(get(handles.motor_m_edit,'String'));
    dm_in = str2num(get(handles.motor_dm_edit,'String'));
    hm_in = str2num(get(handles.motor_h_edit,'String'));
    rm_in = str2num(get(handles.motor_r_edit,'String')); 
    mm = mm_oz/35.274;
    dm = dm_in/39.3701;
    hm = hm_in/39.3701;
    rm = rm_in/39.3701;
    Jx1 = (mm) * ((rm)^2) + (4/3)*(mm) * ((hm)^2) + 2*(mm)*((dm)^2);
    Jy1 = (mm) * ((rm)^2) + (4/3)*(mm) * ((hm)^2) + 2*(mm)*((dm)^2);
    Jz1 = 2*(mm)*((rm)^2) + (4)*(mm)*((dm)^2);
    ms_oz = str2num(get(handles.ESC_m_edit,'String'));
    as_in = str2num(get(handles.ESC_a_edit,'String')); 
    bs_in = str2num(get(handles.ESC_b_edit,'String')); 
    ds_in = str2num(get(handles.ESC_ds_edit,'String')); 
    ms = ms_oz/35.274;
    as = as_in/39.3701;
    bs = bs_in/39.3701;
    ds = ds_in/39.3701;
    Jx2 = ((ms*(as^2))/6) + ((ms*(bs^2))/6) + (2*(ms)*(ds^2));
    Jy2 = ((ms*(as^2))/6) + ((ms*(bs^2))/6) + (2*(ms)*(ds^2));
    Jz2 = ((ms*(as^2 + bs^2))/3) + 4*(ms)*(ds^2);
    mh_oz = str2num(get(handles.HUB_m_edit,'String')); 
    rh_in = str2num(get(handles.HUB_r_edit,'String')); 
    Hh_in = str2num(get(handles.HUB_H_edit,'String')); 
    mh = mh_oz/35.274;
    rh = rh_in/39.3701;
    Hh = Hh_in/39.3701;
    Jx3 = (1/4)*((mh))*((rh)^2) + (1/12)*((mh))*((Hh)^2);
    Jy3 = (1/4)*((mh))*((rh)^2) + (1/12)*((mh))*((Hh)^2);
    Jz3 = (1/2)*((mh))*((rh)^2);
    ma_oz = str2num(get(handles.arms_m_edit,'String'));
    ra_in = str2num(get(handles.arms_r_edit,'String'));
    La_in = str2num(get(handles.arms_L_edit,'String'));
    da_in = str2num(get(handles.arms_da_edit,'String'));
    ma = ma_oz/35.274;
    ra = ra_in/39.3701;
    La = La_in/39.3701;
    da = da_in/39.3701;
    Jx4 = ((3/2)*(ma)*(ra^2)) + ((2/3)*(ma)*(La^2)) + (2)*(ma)*(da^2);
    Jy4 = ((3/2)*(ma)*(ra^2)) + ((2/3)*(ma)*(La^2)) + (2)*(ma)*(da^2);
    Jz4 = ((ma)*((ra)^2)) + (4/3)*((ma) * ((La)^2)) + (4)*(ma)*(da^2);
    totalX = Jx1 + Jx2 + Jx3 + Jx4;
    totalY = Jy1 + Jy2 + Jy3 + Jy4;
    totalZ = Jz1 + Jz2 + Jz3 + Jz4;
    a = num2str(totalX);
    b = num2str(totalY);
    c = num2str(totalZ);
    set(handles.jx,'String',a);
    set(handles.jy,'String',b);
    set(handles.jz,'String',c);
    massTotal = 4*mm + 4*ms + mh + 4*ma;
    set(handles.massTotal_edit,'String',massTotal);
end
guidata(hObject, handles);
function edit39_Callback(hObject, eventdata, handles)
function edit39_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit40_Callback(hObject, eventdata, handles)
function edit40_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
uiload; 
if exist('quadModel')
set(handles.si_button, 'Value', 1)
set(handles.eng_button, 'Value', 0)
set(handles.motor_m_edit,'String',quadModel.motor_m*1000);
set(handles.motor_dm_edit,'String',quadModel.motor_dm*100);
set(handles.motor_h_edit,'String',quadModel.motor_h*100);
set(handles.motor_r_edit,'String',quadModel.motor_r*100);
set(handles.ESC_m_edit,'String',quadModel.ESC_m*1000);
set(handles.ESC_a_edit,'String',quadModel.ESC_a*100);
set(handles.ESC_b_edit,'String',quadModel.ESC_b*100);
set(handles.ESC_ds_edit,'String',quadModel.ESC_ds*100);
set(handles.HUB_m_edit,'String',quadModel.HUB_m*1000);
set(handles.HUB_r_edit,'String',quadModel.HUB_r*100);
set(handles.HUB_H_edit,'String',quadModel.HUB_H*100);
set(handles.arms_m_edit,'String',quadModel.arms_m*1000);
set(handles.arms_r_edit,'String',quadModel.arms_r*100);
set(handles.arms_L_edit,'String',quadModel.arms_L*100);
set(handles.arms_da_edit,'String',quadModel.arms_da*100);
set(handles.ct_edit,'String',quadModel.ct);
set(handles.cq_edit,'String',quadModel.cq);
set(handles.tConstant,'String',quadModel.T);
set(handles.minThr,'String',quadModel.minThr);
set(handles.cr_edit,'String',quadModel.cr);
set(handles.b_edit,'String',quadModel.b);
set(handles.motor_m_unit, 'String', 'g');
set(handles.ESC_m_unit, 'String', 'g');
set(handles.HUB_m_unit, 'String', 'g');
set(handles.Arms_m_unit, 'String', 'g');
set(handles.motor_dm_unit, 'String', 'cm');
set(handles.motor_h_unit, 'String', 'cm');
set(handles.motor_r_unit, 'String', 'cm');
set(handles.ESC_a_unit, 'String', 'cm');
set(handles.ESC_b_unit, 'String', 'cm');
set(handles.ESC_ds_unit, 'String', 'cm');
set(handles.HUB_r_unit, 'String', 'cm');
set(handles.HUB_H_unit, 'String', 'cm');
set(handles.Arms_r_unit, 'String', 'cm');
set(handles.Arms_L_unit, 'String', 'cm');
set(handles.Arms_da_unit, 'String', 'cm');
else
end
function tConstant_Callback(hObject, eventdata, handles)
function tConstant_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes when selected object is changed in unitPanel.
function unitPanel_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue, 'Tag')
    case 'si_button'
        choice1 = questdlg('Changing unit systems erases all cells. Continue?','WARNING','Yes','No','No');
        switch choice1
            case 'Yes'
                set(handles.motor_m_unit, 'String', 'g');
                set(handles.ESC_m_unit, 'String', 'g');
                set(handles.HUB_m_unit, 'String', 'g');
                set(handles.Arms_m_unit, 'String', 'g');
                set(handles.motor_dm_unit, 'String', 'cm');
                set(handles.motor_h_unit, 'String', 'cm');
                set(handles.motor_r_unit, 'String', 'cm');
                set(handles.ESC_a_unit, 'String', 'cm');
                set(handles.ESC_b_unit, 'String', 'cm');
                set(handles.ESC_ds_unit, 'String', 'cm');
                set(handles.HUB_r_unit, 'String', 'cm');
                set(handles.HUB_H_unit, 'String', 'cm');
                set(handles.Arms_r_unit, 'String', 'cm');
                set(handles.Arms_L_unit, 'String', 'cm');
                set(handles.Arms_da_unit, 'String', 'cm');
                set(handles.jx,'String','');
                set(handles.jy,'String','');
                set(handles.jz,'String','');
                set(handles.massTotal_edit,'String','');
                set(handles.ct_edit,'String','');
                set(handles.cq_edit,'String','');
                set(handles.motor_m_edit,'String','');
                set(handles.motor_dm_edit,'String','');
                set(handles.motor_h_edit,'String','');
                set(handles.motor_r_edit,'String','');
                set(handles.ESC_m_edit,'String','');
                set(handles.ESC_a_edit,'String','');
                set(handles.ESC_b_edit,'String','');
                set(handles.ESC_ds_edit,'String','');
                set(handles.HUB_m_edit,'String','');
                set(handles.HUB_r_edit,'String','');
                set(handles.HUB_H_edit,'String','');
                set(handles.arms_m_edit,'String','');
                set(handles.arms_r_edit,'String','');
                set(handles.arms_L_edit,'String','');
                set(handles.arms_da_edit,'String','');
                set(handles.tConstant,'String','');
                set(handles.minThr,'String','');
                set(handles.cr_edit,'String','');
                set(handles.b_edit,'String','');
            case 'No'
                set(handles.eng_button, 'Value', 1)
                set(handles.si_button, 'Value', 0)
                set(handles.motor_m_unit, 'String', 'oz');
                set(handles.ESC_m_unit, 'String', 'oz');
                set(handles.HUB_m_unit, 'String', 'oz');
                set(handles.Arms_m_unit, 'String', 'oz');
                set(handles.motor_dm_unit, 'String', 'in');
                set(handles.motor_h_unit, 'String', 'in');
                set(handles.motor_r_unit, 'String', 'in');
                set(handles.ESC_a_unit, 'String', 'in');
                set(handles.ESC_b_unit, 'String', 'in');
                set(handles.ESC_ds_unit, 'String', 'in');
                set(handles.HUB_r_unit, 'String', 'in');
                set(handles.HUB_H_unit, 'String', 'in');
                set(handles.Arms_r_unit, 'String', 'in');
                set(handles.Arms_L_unit, 'String', 'in');
                set(handles.Arms_da_unit, 'String', 'in');
        end
    case 'eng_button'
        choice2 = questdlg('Changing unit systems erases all cells. Continue?','WARNING','Yes','No','No');
        switch choice2
            case 'Yes'
                set(handles.motor_m_unit, 'String', 'oz');
                set(handles.ESC_m_unit, 'String', 'oz');
                set(handles.HUB_m_unit, 'String', 'oz');
                set(handles.Arms_m_unit, 'String', 'oz');
                set(handles.motor_dm_unit, 'String', 'in');
                set(handles.motor_h_unit, 'String', 'in');
                set(handles.motor_r_unit, 'String', 'in');
                set(handles.ESC_a_unit, 'String', 'in');
                set(handles.ESC_b_unit, 'String', 'in');
                set(handles.ESC_ds_unit, 'String', 'in');
                set(handles.HUB_r_unit, 'String', 'in');
                set(handles.HUB_H_unit, 'String', 'in');
                set(handles.Arms_r_unit, 'String', 'in');
                set(handles.Arms_L_unit, 'String', 'in');
                set(handles.Arms_da_unit, 'String', 'in');
                set(handles.jx,'String','');
                set(handles.jy,'String','');
                set(handles.jz,'String','');
                set(handles.massTotal_edit,'String','');
                set(handles.ct_edit,'String','');
                set(handles.cq_edit,'String','');
                set(handles.motor_m_edit,'String','');
                set(handles.motor_dm_edit,'String','');
                set(handles.motor_h_edit,'String','');
                set(handles.motor_r_edit,'String','');
                set(handles.ESC_m_edit,'String','');
                set(handles.ESC_a_edit,'String','');
                set(handles.ESC_b_edit,'String','');
                set(handles.ESC_ds_edit,'String','');
                set(handles.HUB_m_edit,'String','');
                set(handles.HUB_r_edit,'String','');
                set(handles.HUB_H_edit,'String','');
                set(handles.arms_m_edit,'String','');
                set(handles.arms_r_edit,'String','');
                set(handles.arms_L_edit,'String','');
                set(handles.arms_da_edit,'String','');
                set(handles.tConstant,'String','');
                set(handles.minThr,'String','');
                set(handles.cr_edit,'String','');
                set(handles.b_edit,'String','');
            case 'No'
                set(handles.si_button, 'Value', 1)
                set(handles.eng_button, 'Value', 0)
                set(handles.motor_m_unit, 'String', 'g');
                set(handles.ESC_m_unit, 'String', 'g');
                set(handles.HUB_m_unit, 'String', 'g');
                set(handles.Arms_m_unit, 'String', 'g');
                set(handles.motor_dm_unit, 'String', 'cm');
                set(handles.motor_h_unit, 'String', 'cm');
                set(handles.motor_r_unit, 'String', 'cm');
                set(handles.ESC_a_unit, 'String', 'cm');
                set(handles.ESC_b_unit, 'String', 'cm');
                set(handles.ESC_ds_unit, 'String', 'cm');
                set(handles.HUB_r_unit, 'String', 'cm');
                set(handles.HUB_H_unit, 'String', 'cm');
                set(handles.Arms_r_unit, 'String', 'cm');
                set(handles.Arms_L_unit, 'String', 'cm');
                set(handles.Arms_da_unit, 'String', 'cm');
        end
end
function minThr_Callback(hObject, eventdata, handles)
function minThr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cr_edit_Callback(hObject, eventdata, handles)
function cr_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function b_edit_Callback(hObject, eventdata, handles)
function b_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function arms_da_edit_Callback(hObject, eventdata, handles)
function arms_da_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
