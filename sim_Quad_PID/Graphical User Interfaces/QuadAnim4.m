function varargout = QuadAnim4(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QuadAnim4_OpeningFcn, ...
                   'gui_OutputFcn',  @QuadAnim4_OutputFcn, ...
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
function QuadAnim4_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
handles.AZval = get(handles.AZslider,'Value')*360;
handles.ELval = get(handles.ELslider,'Value')*90;
handles.j     = 1;
handles.skipFlag = 0;
movegui('center')
% --- نمایش لوگو ---
axes(handles.logo_axes); 
try
    imshow('your_logo.jpg'); 
    axis off; 
catch
    axis off; 
end
% ------------------
A = evalin('base', 'yout');
r = .5; d = 1.25; h = .25; % units in inches
a = 1; b = 1; c = 0.2; % units in inches
N = [d  0 h].';
E = [0 -d h].';
W = [0  d h].';
S = [-d 0 h].';
Nr = circlePoints(N, r, 10); Nr = [Nr Nr(:,1)];
Er = circlePoints(E, r, 10); Er = [Er Er(:,1)];
Wr = circlePoints(W, r, 10); Wr = [Wr Wr(:,1)];
Sr = circlePoints(S, r, 10); Sr = [Sr Sr(:,1)];
mN = [d,d;
      0,0;
      h,0];
mE = [0,0;
     -d,-d;
      h,0];
mW = [0,0;
      d,d;
      h,0];
mS = [-d,-d;
       0,0;
       h,0];
bNS = [ d, -d;
        0,  0;
        0,  0];
bEW = [ 0,  0;
        d, -d;
        0,  0];
Top = [ a/2,   0,-a/2,   0;
          0, b/2,   0,-b/2;
        c/2, c/2, c/2, c/2];
Bot = vertcat(Top(1:2,:),-Top(3,:));
NEB = [ a/2, a/2,   0,   0;
          0,   0, b/2, b/2;
        c/2,-c/2,-c/2, c/2];
NWB = [ a/2, a/2,   0,   0;
          0,   0,-b/2,-b/2;
        c/2,-c/2,-c/2, c/2];
SEB = -NWB;
SWB = -NEB;
phi = A(1,4);
the = A(1,5);
psi = A(1,6);
R = [cos(psi)*cos(the) cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi) cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);
       sin(psi)*cos(the) sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi);
       -sin(the)         cos(the)*sin(phi)                            cos(the)*cos(phi)];
U = A(:,7);
V = A(:,8);
W = A(:,9);
Vi = zeros(length(A),3);
MvMax= max(sqrt(U.^2+V.^2+W.^2));
Vb = 3/MvMax*[U(1), V(1), W(1)]';
Vi(1,:) = R*Vb;
if (~evalin('base','quadModel.plusConfig'))
    Rz = [ sqrt(2)/2, sqrt(2)/2, 0;
          -sqrt(2)/2,sqrt(2)/2, 0;
                   0,          0, 1];
    Nr = Rz*Nr;
    Er = Rz*Er;
    Wr = Rz*Wr;
    Sr = Rz*Sr;
    mN = Rz*mN;
    mE = Rz*mE;
    mW = Rz*mW;
    mS = Rz*mS;
    bNS = Rz*bNS;
    bEW = Rz*bEW;
    Top = Rz*Top;
    Bot = Rz*Bot;
    NEB = Rz*NEB;
    NWB = Rz*NWB;
    SWB = Rz*SWB;
    SEB = Rz*SEB;
end
NrR = R*Nr;
ErR = R*Er;
WrR = R*Wr;
SrR = R*Sr;
mNr = R*mN;
mEr = R*mE;
mWr = R*mW;
mSr = R*mS;
bNSR = R*bNS;
bEWR = R*bEW;
TopR = R*Top;
BotR = R*Bot;
NEBR = R*NEB;
NWBR = R*NWB;
SWBR = R*SWB;
SEBR = R*SEB;
axes(handles.axes1)
plot3(bNSR(1,:),bNSR(2,:),bNSR(3,:),'b','LineWidth',3)
hold on
plot3(bEWR(1,:),bEWR(2,:),bEWR(3,:),'b','LineWidth',3)
plot3(mNr(1,:),mNr(2,:),mNr(3,:),'k','LineWidth',4)
plot3(mEr(1,:),mEr(2,:),mEr(3,:),'k','LineWidth',4)
plot3(mWr(1,:),mWr(2,:),mWr(3,:),'k','LineWidth',4)
plot3(mSr(1,:),mSr(2,:),mSr(3,:),'k','LineWidth',4)
plot3(NrR(1,:),NrR(2,:),NrR(3,:),'g')
plot3(ErR(1,:),ErR(2,:),ErR(3,:),'k')
plot3(WrR(1,:),WrR(2,:),WrR(3,:),'k')
plot3(SrR(1,:),SrR(2,:),SrR(3,:),'g')
grey = [0.5 0.5 0.5];
top = fill3(TopR(1,:),TopR(2,:),TopR(3,:),'r'); alpha(top,0.8);
bot = fill3(BotR(1,:),BotR(2,:),BotR(3,:),'g'); alpha(bot,0.8);
ne  = fill3(NEBR(1,:),NEBR(2,:),NEBR(3,:),'c'); alpha(ne,0.8);
nw  = fill3(NWBR(1,:),NWBR(2,:),NWBR(3,:),grey); alpha(nw,0.8);
sw  = fill3(SWBR(1,:),SWBR(2,:),SWBR(3,:),grey); alpha(sw,0.8);
se  = fill3(SEBR(1,:),SEBR(2,:),SEBR(3,:),grey); alpha(se,0.8);
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
view(handles.AZval,handles.ELval)
grid on
omi = zeros(length(A),3);
P = A(:,1); Q = A(:,2); Rw = A(:,3);
MombMax = max(sqrt(P.^2+Q.^2+Rw.^2));
omb = 3/MombMax*[P,Q,Rw].';
omi(1,:) = R*omb(:,1);
qp1 = quiver3(0,0,0,omi(1,1),omi(1,2),omi(1,3),'ro');
qp2 = quiver3(0,0,0,Vi(1,1),Vi(1,2),Vi(1,3),'k');
hold off
minX = min(A(:,10));
minY = min(A(:,11));
maxX = max(A(:,10));
maxY = max(A(:,11));
maxZ = max(A(:,12));
axes(handles.axes2)
X = A(1:1,10); Y = A(1:1,11); Z = A(1:1,12);
scatter3(X,Y,Z,36,'blue')
hold on
fill3([minX-1 maxX+1 maxX+1 minX-1],...
    [minY-1 minY-1 maxY+1 maxY+1],...
    [0 0 0 0],'g');
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
xlim([minX-1 maxX+1])
ylim([minY-1 maxY+1])
zlim([-0.1 maxZ+1])
view(handles.AZval,handles.ELval)
axis square
grid on
hold off
handles.output = hObject;
guidata(hObject, handles);
function varargout = QuadAnim4_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
function startStop_Callback(hObject, eventdata, handles)
plotAnim(hObject, handles)
function ELslider_Callback(hObject, eventdata, handles)
handles.ELval = get(hObject,'Value')*90;
handles.AZval = get(handles.AZslider,'Value')*360;
view(handles.axes1,handles.AZval,handles.ELval);
view(handles.axes2,handles.AZval,handles.ELval);
drawnow
guidata(hObject,handles)
function ELslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function AZslider_Callback(hObject, eventdata, handles)
handles.AZval = get(hObject,'Value')*360;
handles.ELval = get(handles.ELslider,'Value')*90;
view(handles.axes1,handles.AZval,handles.ELval);
view(handles.axes2,handles.AZval,handles.ELval);
drawnow
guidata(hObject,handles);
function AZslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function axes1_CreateFcn(hObject, eventdata, handles)
function points = circlePoints(center, radius, numberOfPoints)
c = center.';
r = radius;
n = numberOfPoints;
th = (0:n-1)'/n*2*pi;
x = r*cos(th) + c(1);
y = r*sin(th) + c(2);
points = [x,y].';
    if length(c) > 2
        z = ones(size(x))*c(3);
        points = [x, y, z].';
    end
function frameSkips_Callback(hObject, eventdata, handles)
handles.frameSkipVal = str2double(get(hObject,'String'));
guidata(hObject,handles);
function frameSkips_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.frameSkipVal = str2double(get(hObject,'String'));
guidata(hObject,handles);
function plotAnim(button, handles)
hObject = button;
mode = get(hObject,'String');
handles = guidata(gcbo);
j = handles.j;
if strcmp(mode,'Start')
    set(hObject,'String','Stop')
    guidata(hObject,handles)
    A = evalin('base', 'yout');
    tout = evalin('base', 'tout');
    frameSkipVal = str2double(get(handles.frameSkips,'String'))+1;
    
    r = .5; d = 1.25; h = .25; % units in inches
    a = 1; b = 1; c = 0.2; % units in inches
    N = [d  0 h].';
    E = [0 -d h].';
    W = [0  d h].';
    S = [-d 0 h].';
    Nr = circlePoints(N, r, 20); Nr = [Nr Nr(:,1)];
    Er = circlePoints(E, r, 20); Er = [Er Er(:,1)];
    Wr = circlePoints(W, r, 20); Wr = [Wr Wr(:,1)];
    Sr = circlePoints(S, r, 20); Sr = [Sr Sr(:,1)];
    mN = [d,d;
          0,0;
          h,0];
    mE = [0,0;
         -d,-d;
          h,0];
    mW = [0,0;
          d,d;
          h,0];
    mS = [-d,-d;
           0,0;
           h,0];
    bNS = [ d, -d;
            0,  0;
            0,  0];
    bEW = [ 0,  0;
            d, -d;
            0,  0];
    Top = [ a/2,   0,-a/2,   0;
              0, b/2,   0,-b/2;
            c/2, c/2, c/2, c/2];
    Bot = vertcat(Top(1:2,:),-Top(3,:));
    NEB = [ a/2, a/2,   0,   0;
          0,   0, b/2, b/2;
        c/2,-c/2,-c/2, c/2];
    NWB = [ a/2, a/2,   0,   0;
              0,   0,-b/2,-b/2;
            c/2,-c/2,-c/2, c/2];
    SEB = -NWB;
    SWB = -NEB;
    phi = A(1,4);
    the = A(1,5);
    psi = A(1,6);
    R = [cos(psi)*cos(the) cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi) cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);
           sin(psi)*cos(the) sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi);
           -sin(the)         cos(the)*sin(phi)                            cos(the)*cos(phi)];
    
    if (~evalin('base','quadModel.plusConfig'))
        Rz = [ sqrt(2)/2, sqrt(2)/2, 0;
               -sqrt(2)/2,sqrt(2)/2, 0;
                       0,          0, 1];
        Nr=Rz*Nr;
        Er=Rz*Er;
        Wr=Rz*Wr;
        Sr=Rz*Sr;
        mN = Rz*mN;
        mE = Rz*mE;
        mW = Rz*mW;
        mS = Rz*mS;
        bNS = Rz*bNS;
        bEW = Rz*bEW;
        Top = Rz*Top;
        Bot = Rz*Bot;
        NEB = Rz*NEB;
        NWB = Rz*NWB;
        SWB = Rz*SWB;
        SEB = Rz*SEB;
    end
    
    U = A(:,7);
    V = A(:,8);
    W = A(:,9);
    Vi = zeros(length(A),3);
    MvMax= max(sqrt(U.^2+V.^2+W.^2));
    Vb = 3/MvMax*[U, V, W].';
    Vi(1,:) = R*Vb(:,1);
    P = A(:,1);
    Q = A(:,2);
    Rw = A(:,3);
    omi = zeros(length(A),3);
    MombMax = max(sqrt(P.^2+Q.^2+Rw.^2));
    omb = 3/MombMax*[P,Q,Rw].';
    omi(1,:) = R*omb(:,1);
    phi = A(:,4);
    the= A(:,5);
    psi = A(:,6);
    minX = min(A(:,10));
    minY = min(A(:,11));
    maxX = max(A(:,10));
    maxY = max(A(:,11));
    maxZ = max(A(:,12));
    
    R = cell(length(A),1);
    for i = 1:length(A)
    R{i,1} = [cos(psi(i))*cos(the(i)) cos(psi(i))*sin(the(i))*sin(phi(i))-sin(psi(i))*cos(phi(i)) cos(psi(i))*sin(the(i))*cos(phi(i))+sin(psi(i))*sin(phi(i));
              sin(psi(i))*cos(the(i)) sin(psi(i))*sin(the(i))*sin(phi(i))+cos(psi(i))*cos(phi(i)) sin(psi(i))*sin(the(i))*cos(phi(i))-cos(psi(i))*sin(phi(i));
              -sin(the(i))         cos(the(i))*sin(phi(i))                            cos(the(i))*cos(phi(i))];
    end
    while ( strcmp(get(hObject,'String'),'Stop'))
        handles = guidata(gcbo);
        j = handles.j;
        
        guidata(hObject,handles);
        Vi(j,:) = R{j,1}*Vb(:,j);
        NrR = R{j,1}*Nr;
        ErR = R{j,1}*Er;
        WrR = R{j,1}*Wr;
        SrR = R{j,1}*Sr;
        bNSR = R{j,1}*bNS;
        bEWR = R{j,1}*bEW;
        TopR = R{j,1}*Top;
        BotR = R{j,1}*Bot;
        mNr = R{j,1}*mN;
        mEr = R{j,1}*mE;
        mWr = R{j,1}*mW;
        mSr = R{j,1}*mS;
        NEBR = R{j,1}*NEB;
        NWBR = R{j,1}*NWB;
        SWBR = R{j,1}*SWB;
        SEBR = R{j,1}*SEB;
        axes(handles.axes1)
        plot3(bNSR(1,:),bNSR(2,:),bNSR(3,:),'b','LineWidth',3)
        hold on
        plot3(bEWR(1,:),bEWR(2,:),bEWR(3,:),'b','LineWidth',3)
        plot3(NrR(1,:),NrR(2,:),NrR(3,:),'g')
        plot3(ErR(1,:),ErR(2,:),ErR(3,:),'k')
        plot3(WrR(1,:),WrR(2,:),WrR(3,:),'k')
        plot3(SrR(1,:),SrR(2,:),SrR(3,:),'g')
        plot3(mNr(1,:),mNr(2,:),mNr(3,:),'k','LineWidth',4)
        plot3(mEr(1,:),mEr(2,:),mEr(3,:),'k','LineWidth',4)
        plot3(mWr(1,:),mWr(2,:),mWr(3,:),'k','LineWidth',4)
        plot3(mSr(1,:),mSr(2,:),mSr(3,:),'k','LineWidth',4)
        fill3(TopR(1,:),TopR(2,:),TopR(3,:),'r'); alpha(0.8);
        fill3(BotR(1,:),BotR(2,:),BotR(3,:),'g'); alpha(0.8);
        grey = [0.5 0.5 0.5];
        ne  = fill3(NEBR(1,:),NEBR(2,:),NEBR(3,:),'c'); alpha(ne,0.8);
        nw  = fill3(NWBR(1,:),NWBR(2,:),NWBR(3,:),grey); alpha(nw,0.8);
        sw  = fill3(SWBR(1,:),SWBR(2,:),SWBR(3,:),grey); alpha(sw,0.8);
        se  = fill3(SEBR(1,:),SEBR(2,:),SEBR(3,:),grey); alpha(se,0.8);
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        xlim([-2 2])
        ylim([-2 2])
        zlim([-2 2])
        handles.AZval = get(handles.AZslider,'Value')*360;
        handles.ELval = get(handles.ELslider,'Value')*90;
        view(handles.axes1,handles.AZval,handles.ELval)
        omi(j,:) = R{j,1}*omb(:,j);
        qp1 = quiver3(0,0,0,omi(j,1),omi(j,2),omi(j,3),'r');
        qp2 = quiver3(0,0,0,Vi(j,1),Vi(j,2),Vi(j,3),'k');
        axis square
        grid on
        hold off
        drawnow
        
        if (j==1)
            cla(handles.axes2)
        else if (handles.skipFlag==1)
                cla(handles.axes2)
                X = A(1:frameSkipVal:j,10);
                Y = A(1:frameSkipVal:j,11);
                Z = A(1:frameSkipVal:j,12);
                axes(handles.axes2)
                hold on
                scatter3(X,Y,Z,36,'blue');
            end
        end
        X = A(j,10); Y = A(j,11); Z = A(j,12);
        axes(handles.axes2)
        hold on
        scatter3(X,Y,Z,36,'blue');
        if (j == 1 || handles.skipFlag==1)
            fill3([minX-1 maxX+1 maxX+1 minX-1],...
                  [minY-1 minY-1 maxY+1 maxY+1],...
                  [0 0 0 0],'g');
            alpha(0.7);
            handles.skipFlag = 0;
            guidata(hObject,handles)
        end
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        xlim([minX-1 maxX+1])
        ylim([minY-1 maxY+1])
        zlim([-0.1 maxZ+1])
        handles.AZval = get(handles.AZslider,'Value')*360;
        handles.ELval = get(handles.ELslider,'Value')*90;
        view(handles.axes2,handles.AZval,handles.ELval)
        axis square
        grid on
        drawnow
        hold off
        if (j == size(A,1))
            handles.j = 1;
            set(hObject,'String','Start');
            guidata(hObject,handles);
            return
            
        else if (j+frameSkipVal<size(A,1))
        j = j+frameSkipVal;
            else
            j = size(A,1);
            end
        end
        set(handles.simTime,'String',num2str(tout(j)));
        handles = guidata(gcbo);
        if (handles.skipFlag == 1)
            guidata(hObject,handles)
            return
        else
        handles.j = j;
        guidata(hObject,handles);
        end
    end
    else if strcmp(get(hObject,'String'),'Stop')
        set(hObject,'String','Start');
        end
end
function simTime_CreateFcn(hObject, eventdata, handles)
function axes2_CreateFcn(hObject, eventdata, handles)
function XYview_Callback(hObject, eventdata, handles)
set(handles.ELslider,'Value',1);
set(handles.AZslider,'Value',0);
view(handles.axes1,0,90);
view(handles.axes2,0,90);
drawnow
guidata(hObject,handles);
function XZview_Callback(hObject, eventdata, handles)
set(handles.ELslider,'Value',0);
set(handles.AZslider,'Value',0);
view(handles.axes1,0,0);
view(handles.axes2,0,0);
drawnow
guidata(hObject,handles);
function YZview_Callback(hObject, eventdata, handles)
set(handles.ELslider,'Value',0);
set(handles.AZslider,'Value',90/360);
view(handles.axes1,90,0);
view(handles.axes2,90,0);
drawnow
guidata(hObject,handles);
function defaultView_Callback(hObject, eventdata, handles)
set(handles.ELslider,'Value',25/90);
set(handles.AZslider,'Value',35/360)
view(handles.axes1,35,25);
view(handles.axes2,35,25);
drawnow
guidata(hObject,handles);
function timeSkipButton_Callback(hObject, eventdata, handles)
set(handles.startStop,'String','Start');
tout = evalin('base', 'tout');
TimeSkip = str2double(get(handles.timeSkipEditBox,'String'));
handles.j = find(tout<=TimeSkip, 1,'last');
handles.skipFlag = 1;
guidata(hObject,handles);
function timeSkipEditBox_Callback(hObject, eventdata, handles)
handles.TimeSkip = str2double(get(hObject,'String'));
handles.skipFlag = 0;
guidata(hObject,handles);
function timeSkipEditBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.TimeSkip = 0;
guidata(hObject,handles);
function figure1_CreateFcn(hObject, eventdata, handles)
