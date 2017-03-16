function varargout = BCOMS(varargin)
% BCOMS MATLAB code for BCOMS.fig
%      BCOMS, by itself, creates a new BCOMS or raises the existing
%      singleton*.
%
%      H = BCOMS returns the handle to a new BCOMS or the handle to
%      the existing singleton*.
%
%      BCOMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BCOMS.M with the given input arguments.
%
%      BCOMS('Property','Value',...) creates a new BCOMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BCOMS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BCOMS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BCOMS

% Last Modified by GUIDE v2.5 13-Feb-2017 11:13:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BCOMS_OpeningFcn, ...
                   'gui_OutputFcn',  @BCOMS_OutputFcn, ...
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


% --- Executes just before BCOMS is made visible.
function BCOMS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BCOMS (see VARARGIN)

% Initialize
handles.flagRoiExtract = 0;
handles.flagEmbReg = 0;
handles.flagMemb = 0;
handles.memraneDir = './';

% Choose default command line output for BCOMS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BCOMS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BCOMS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% Image data
% pushbutton
function pushbutton5_Callback(hObject, eventdata, handles)
[fileName, pathName] = uigetfile('*','Select membrane image file');
handles.memraneDir = pathName;
handles.memraneFilename = [pathName, fileName];
if handles.memraneFilename == 0; return; end
set(handles.edit1, 'String', num2str(handles.memraneFilename));
guidata(hObject, handles);

function pushbutton6_Callback(hObject, eventdata, handles)
[fileName, pathName] = uigetfile('*','Select nuclear segmentation file', handles.memraneDir);
handles.nucleusFilename = [pathName, fileName];
if handles.nucleusFilename == 0; return; end
set(handles.edit2, 'String', num2str(handles.nucleusFilename));
guidata(hObject, handles);

function pushbutton7_Callback(hObject, eventdata, handles)
handles.resultDir = uigetdir(handles.memraneDir,'Select results directory');
if handles.resultDir == 0; return; end
set(handles.edit3, 'String', num2str(handles.resultDir));
guidata(hObject, handles);

function pushbutton1_Callback(hObject, eventdata, handles)


% edit
function edit1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

handles.memraneFilename = get(hObject,'String');
guidata(hObject, handles);

function edit1_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
handles.nucleusFilename = get(hObject,'String');
guidata(hObject, handles);

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)
handles.resultDir = get(hObject,'String');
guidata(hObject, handles);

function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Image information

function pushbutton2_Callback(hObject, eventdata, handles)

function edit4_Callback(hObject, eventdata, handles)
handles.resXY = str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)
handles.resZ = str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
handles.numZ = str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)
handles.numT = str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% ROI extraction
function pushbutton8_Callback(hObject, eventdata, handles)

if ~isfield(handles, 'resXY') || ~isfield(handles, 'resZ') || ~isfield(handles, 'numZ') || ~isfield(handles, 'numT')
    ed = errordlg('Please set the Image information in advance','Error');
    set(ed, 'WindowStyle', 'modal');
    uiwait(ed);
    return
end

% handles.ROIDir = [handles.resultDir, '\ROI'];

% roiSelect( 'MyCele', handles.memraneFilename, [], handles.ROIDir, 1, handles.numT, handles.numZ, 1 );

% h = msgbox('Extracting ROI region');
% extractImg( 'MyCele', handles.memraneFilename, handles.nucleusFilename, 1, handles.numT, handles.numZ, handles.resultDir, 1 );
% delete(h);

h = msgbox('Reading image files');
handles.membImgDir = [handles.resultDir, '\MembraneImage'];
handles.nucImgDir = [handles.resultDir, '\NuclearImage'];

tifRead(handles.memraneFilename, handles.numZ, handles.numT, handles.membImgDir)
tifRead(handles.nucleusFilename, handles.numZ, handles.numT, handles.nucImgDir)

handles.flagRoiExtract = 1;
guidata(hObject, handles);
delete(h);
h = msgbox('Finished reading');

%% Embryonic region segmentation

function pushbutton3_Callback(hObject, eventdata, handles)
if handles.flagRoiExtract == 0
    ed = errordlg('Please Run the ROI extraction in advance','Error');
    set(ed, 'WindowStyle', 'modal');
    uiwait(ed);
    return
end

handles.membImgROIDir = [handles.resultDir, '\Membrane\ImgRoi'];
handles.nucSegROIDir = [handles.resultDir, '\Nucleus\ImgRoi'];
handles.roiExtractedDir = [handles.resultDir, '\ROI\Extracted'];
handles.embRegDir = [handles.resultDir, '\EmbReg'];

h = msgbox('Computing embryonic region');
ME = 0;
try
    embryonicRegion(handles.membImgDir, handles.nucImgDir, handles.embRegDir, handles.volRatioThresh);
catch ME
    delete(h);
    rethrow(ME)
end
delete(h);
if ME==0
    h = msgbox('Embryonic region segmentation completed');
end

handles.flagEmbReg = 1;
guidata(hObject, handles);

function edit8_Callback(hObject, eventdata, handles)
handles.volRatioThresh = str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.volRatioThresh = str2double(get(hObject,'String'));
guidata(hObject, handles);

%% Cell membrane segmentation
function pushbutton4_Callback(hObject, eventdata, handles)
if handles.flagRoiExtract == 0
    ed = errordlg('Please Run the ROI extraction in advance','Error');
    set(ed, 'WindowStyle', 'modal');
    uiwait(ed);
    return
end
if handles.flagEmbReg == 0
    ed = errordlg('Please Run the Embryonic region segmentation in advance','Error');
    set(ed, 'WindowStyle', 'modal');
    uiwait(ed);
    return
end

handles.membSegDir = [handles.resultDir, '\MembraneSegmentation'];
handles.embRegStackDir = [handles.resultDir, '\EmbReg\Stack'];

h = msgbox('Computing membrane segmentation');
ME = 0;
try
    waterMembrane( handles.membImgDir, handles.nucImgDir, handles.embRegStackDir, handles.membSegDir, handles.resXY, handles.resZ );
%     simpleWater( handles.membImgDir, handles.nucImgDir, handles.embRegStackDir, handles.membSegDir, handles.resXY, handles.resZ );
catch ME
    delete(h);
%     e = errordlg(ME.message);
    rethrow(ME)
end
delete(h);
if ME==0
    h = msgbox('Membrane segmentation completed');
end
