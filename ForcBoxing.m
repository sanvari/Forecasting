function varargout = ForcBoxing(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ForcBoxing_OpeningFcn, ...
                   'gui_OutputFcn',  @ForcBoxing_OutputFcn, ...
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


% --- Executes just before ForcBoxing is made visible.
function ForcBoxing_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
global teta
teta = 0.3;
    
% --- Outputs from this function are returned to the command line.
function varargout = ForcBoxing_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;


function [ACF, Lags, Bounds]=SE_ACF (iData)
% This function calculates ACF for input data
% Input:
%     Series: input series
% Output:
%     ACF: ACF for input series
%     Lags: ACF Lags
%     Bounds: ACF bounds
[ACF, Lags, Bounds]=autocorr (iData);
Lags = Lags (2:length(Lags));
ACF = ACF (2:length(ACF));


function [PACF, Lags, bounds]=SE_PACF(iData)
% This function calculates PACF for input data
% Input:
%     Series: input series
% Output:
%     PACF: PACF for input series
%     Lags: PACF Lags
%     bounds: PACF bounds
[partialACF, Lags, bounds]=parcorr(iData);
PACF=partialACF(2:length(partialACF));
Lags=Lags(2:length(Lags));


function [fullPath] = SelectFile()
[fname, fpath] = uigetfile({'*.txt;*.csv','Text files(*.txt, *.csv)'
                     '*.*', 'All Files (*.*)'}, 'Select the input file');
fullPath = fullfile(fpath, fname);

global df
U = load(fullPath);
inputData = U(:,1);
%disp('original serie');
%disp(inputData);
figure(20)
plot(inputData);
[rows, cols] = size(inputData);
xlim([0, rows])
set(gca,'XTick', 0:10:rows)
set(gca,'XTickLabel', 0:10:rows)
title('Input Data');

%X = var(inputData);
%figure(22);
%plot(X);

df = Differecings(inputData);

df.LogTrans();
df.seasonaldiff();
df.NormalDiff();
figure(21);
subplot(2,1,1);
autocorr(df.finalSeries);
title('ACF of Stationary Data')
subplot(2,1,2);
parcorr(df.finalSeries);
title('PACF of Stationary Data');

function txtFilePath_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btnSelectFile.
function btnSelectFile_Callback(~, ~, handles)
fPath = SelectFile();
set(handles.txtFilePath, 'string', fPath);

% --- Executes on button press in AutoCorrelation.
function AutoCorrelation_Callback(~, ~, ~)
global df 
figure(1)
plot(df.oldSeries);

%nnn = df.SeriesSize;
xlim([0, df.SeriesSize])
set(gca,'XTick', 0:10:df.SeriesSize)
set(gca,'XTickLabel', 0:10:df.SeriesSize)
title('Passenger Demand Monday Mornings Merter Station: 1July till 27 November 2012, 7-10 am')

% --- Executes on button press in PAutoCorrelation.
function PAutoCorrelation_Callback(~, ~, ~)
global df 
figure(2);
subplot(2,1,1);
autocorr(df.oldSeries);
subplot(2,1,2);
parcorr(df.oldSeries);

% --- Executes on button press in DAutoCorrelation.
function DAutoCorrelation_Callback(~, ~, ~)
global df 
figure(3);
plot(df.finalSeries);
xlim([0,df.SeriesSize]);
set(gca,'XTick', 0:10:df.SeriesSize);
set(gca,'XTickLabel', 0:10:df.SeriesSize);
title('Differenced Log Quarterly Australian CPI');

% --- Executes on button press in DPAutoCorrelation.
function DPAutoCorrelation_Callback(~, ~, ~)
global df 
figure(4);
subplot(2,1,1);
autocorr(df.finalSeries);
subplot(2,1,2);
parcorr(df.finalSeries);

% --- Executes during object creation, after setting all properties.
function txtP_CreateFcn(hObject, ~, ~)
% hObject    handle to txtP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function txtQ_CreateFcn(hObject, ~, ~)
% hObject    handle to txtQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function txtI_CreateFcn(hObject, ~, ~)
% hObject    handle to txtI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
                                      
% --------------------------------------------------------------------
function menuSelectFile_Callback(~, ~, handles)
%fPath = SelectFile();
%set(handles.txtFilePath, 'string', fPath);
set(handles.edit5, 'string', num2str(223));

% --------------------------------------------------------------------
function menuExit_Callback(~, ~, handles)
close(handles.frmModel);
   

function btnForecast_Callback(hObject, eventdata, handles)
global df
myfrc= MyForecast(df);
myfrc.DoForecast(df);