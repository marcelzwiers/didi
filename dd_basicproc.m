function varargout = dd_basicproc(varargin)

% DD_BASICPROC converts DICOM data to Nifti files, corrects for subject
% and/or cardiac motion artefacts, normalises the data to the T2-EPI
% ICBM-template, and computes the diffusion tensor and its main derivative
% scalar measures.
%
%
% Requirements:
% ~~~~~~~~~~~~~
% - SPM12
% - Diffusion toolbox (http://sourceforge.net/projects/spmtools)
% - BET (http://www.fmrib.ox.ac.uk/fsl/bet)
% - Helper functions in the directory where this file is/was found
% - Distributed computing toolbox (optional). See:
%   http://fieldtrip.fcdonders.nl/development/qsub
%
%
% GUI Usage:
% ~~~~~~~~~~
% Running DD_BASICPROC by itself (no input arguments) creates a new gui
% or raises the existing one. This gui has three main panels:
%
% 1) DATA PANEL (blue). Here the input data can be given as raw DICOM
%    files, or as (previously processed) Nifti files. This latter option is
%    particularly usefull when the automatic procedure (e.g. the brain
%    masking or realignment) failed in individual data sets: just correct
%    these cases manually and use this option to rerun the job selecting
%    the corrected (Nifti) files. For example, when you want to correct the
%    brain mask, just redo the brain extraction manually where needed (e.g.
%    using BET or MRICro) and rerun your job with the Nifti option selected.
%    Another example would be to redo the realignment in some data-sets
%    (e.g. using SPM) and rerun your job with the Nifti-option, selecting
%    the updated 'r*.img' files and setting the realignment-option to
%    'none' in the processing panel (as these files are already realigned).
%       Note that multiple data sets (subjects) can be added to construct
%    large batch jobs. The data-set (nr) selector can be used to inspect
%    and remove your added data sets ('remove' removes the selected data
%    set only).
%
% 2) PROCESSING PANEL (orange). Various options for how the DTI processing
%    should be performed can be selected here.
%		First, the raw data can be denoised using the LPCA filter (recommended)
%	 as described in: J. V. Manjon, P. Coupe, L. Concha, A. Buades, D. L.
%	 Collins, M. Robles (2013). Diffusion Weighted Image Denoising using
%	 overcomplete Local PCA. PLoS ONE 8: doi: 10.1371/journal.pone.0073021. For
%	 other filters see: https://sites.google.com/site/pierrickcoupe/softwares/
%    denoising-for-medical-imaging/dwi-denoising/dwi-denoising-software
%		Using the 'artefact detection' options, data sets can be inspected and
%	 corrected for (cardiac) motion artefacts interactively or automatically.
%       The next processing step concerns correction for geometric
%    distortions ('Eddy-current correction') and subject motion. Various
%    options are possible:
%    - 'sq DT error': This correction methods corrects for both the eddy
%      currents and subject motion simultaneously, and is based on
%      minimization of the residual diffusion tensor error. In practice, this
%	   method is often not that effective (ref: Andersson JLR, Skare S (2002). A
%	   Model-Based Method for Retrospective Correction of Geometric Distortions
%	   in Diffusion-Weighted EPI. NeuroImage 16: 177-199).
%    - 'mutual info': This affine transformation is the most robust option but
%	   computationally much more intensive. This method uses the SPM
%	   coregistration function, and is based on minization of the normalised
%	   mutual information measure (see help('spm_coreg') or
%	   spm_help('jobs.spatial.coreg')).
%    - 'sq error': An option to quickly realign using the SPM realign
%      function (see spm_realign). This rigid body method is based on
%      minimization of difference in image intensity and is the least
%      precise method of the three, and does not correct for eddy-current
%      distortions (to be implemented in future versions).
%    The next processing step is (non-linear) normalization on the SPM
%    (ICBM) T2-EPI template. The normalisation can include the diffusion
%    weighted images ('DW images') or only the computes DTI measures
%    ('results only'). This latter option normalises the FA, RA and MD/ADC
%    measures (see further below), but not the eigenvectors and values - for
%    that you should normalise using 'DW images'. In general, it is
%    probably best to use the 'results only' option when doing a conventional
%    group study on DTI measures. Normalisation is more reliable when the
%    T1 reference scan option is used (see further below).
%       Next, an average b0/DWI-mask can be automatically constructed (using
%    BET) and used to mask your results. The BET fractional intensity threshold
%    can be set to a customized value between 0 and 1 (a smaller value gives a
%    larger brain outline estimate). Always check automatically constructed
%    masks (e.g. with fslview; load the b0-image and then load the mask as an
%    image overlay). For many operations, such as artefact detection or
%    realignment, a mask is created (but not necessarily applied to your
%    results).
%       If you have also collected a T1 reference image you can coregister
%    the diffusion weighted images to this reference image (select 'use T1
%    reference scan'). If you choose to normalise the diffusion weighted
%    images, the (coregistered) reference image will be used to improve
%    normalisation: The spatial normalisation parameters will be computed
%    on the T1 reference image (instead of on the average b0 image) and
%    applied to the diffusion weighted images. Note that this should
%    produce a more reliable normalization.
%		The diffusion tensors can be estimated using the PATCH method (see
%	 patch.m or doi:10.1016/j.neuroimage.2010.06.014) directly or using the
%	 artefact corrected data in the SPM diffusion toolbox.
%       In the final processing steps, the DTI measures are computed. The
%    possible options are:
%    - 'fractional anisotropy (FA)': Most commonly used scalar DTI measure
%      that expresses the diffusion tensor directionality.
%    - 'relative anistropy (RA)': Alternative measure for the FA.
%    - 'mean diffusivity (MD/ADC)': Most used scalar DTI measure for that
%      expresses the size of the diffusion tensor with important clinical
%      applications (e.g. edema/stroke detection).
%    - 'eigenvectors & values': These can be used for various purposes, for
%      example to visualise the local direction of main diffusion (see
%      the SPM Diffusion toolbox).
%    You don't need to select any of these options if you are only
%    interested in doing fiber tractography.
%       If you want to do fiber tractography using the FSL or Camino software
%    package, you can use the export options to create appropriate files (in
%    FDT_Data and in Camino_Data, respectively).
%
% 3) JOB PANEL (green). The buttons in this panel can be used to save and run
%    the current job ('Save & Run') or to load the settings of a previously
%    saved job-file ('Load'; see DATA PANEL above). The results can be viewed
%    using SPM's display or checkreg functionality (the latter has
%    DTI/spm_orthview plug-ins that can be used by right-clicking the image).
% 
% Commandline Usage:
% ~~~~~~~~~~~~~~~~~~
% Previously saved jobs can also be run directly from the Matlab
% commandline or from witin m-files, without gui interaction. If the
% variable Job is not in your workspace (the Job structure is put in the
% caller's workspace when a gui-job is run), or if you want to run a
% different job than the last one, load this variable in your workspace,
% e.g. by running:
%
% >> load('[the name of your job.mat-file]')
%
% Then, simply run:
%
% >> dd_basicproc('jobman', Job, '[the name of your job.mat-file]')
%
% NB: 1) Commandline-usage is supported for Matlab 7.4 (R2007a) and higher only.
%	  2) Make sure your Matlab-path is correct: commandline usage skips the
%	     normal path-checks
%
% _________________________________________________________________________
% DD_BASICPROC was developed by Marcel Zwiers. Please report all bugs to:
% m.zwiers@donders.ru.nl
%
% Ver 4.0, 19/12/2014
%
% See also: DD_PATCH

% DD_BASICPROC M-file for dd_basicproc.fig
%      DD_BASICPROC, by itself, creates a new DD_BASICPROC or raises the existing
%      singleton*.
%
%      H = DD_BASICPROC returns the handle to a new DD_BASICPROC or the handle to
%      the existing singleton*.
%
%      DD_BASICPROC('Property','Value',...) creates a new DD_BASICPROC using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to dd_basicproc_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DD_BASICPROC('CALLBACK') and DD_BASICPROC('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DD_BASICPROC.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Edit the above text to modify the response to help dd_basicproc

% Last Modified by GUIDE v2.5 12-Apr-2019 19:25:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dd_basicproc_OpeningFcn, ...
                   'gui_OutputFcn',  @dd_basicproc_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before dd_basicproc is made visible.
function dd_basicproc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Check if we are all set
% if ~verLessThan('matlab', '7.13.0')
% 	warning('You are using Matlab %s. Only versions less recent than 7.13 are supported (due to obsolete fileparts output arguments)', version)
% end
if ~any(strcmpi(spm('ver'), {'SPM12', 'SPM12b'}))
    Ans = questdlg('This toolbox works with SPM12 or higher only! Do you want to set your path to SPM12?', 'SPM version conflict');
	if strcmp(Ans, 'Yes')
		pathtool
	else
		return
	end
end
WhatdTB = what('Diffusion');  % Check if the Diffusion-toolbox on the Matlab-path
xTB     = spm('TBs');         % Check if the Diffusion-toolbox is installed in SPM
dTB     = find(strcmp('Diffusion', {xTB.name}));
if numel(WhatdTB)==0
    if isempty(dTB), error('The SPM Diffusion-toolbox is not installed'), end
	disp('=> Adding the Diffusion toolbox to your path')
    addpath(genpath(xTB(dTB).dir))
elseif numel(WhatdTB)==1 && (isempty(dTB) || ~strcmp(xTB(dTB).dir, WhatdTB.path))
    disp(['Unexpected Diffusion-toolbox path: ' WhatdTB.path])
elseif numel(WhatdTB)>1
    warndlg(sprintf(['More than 1 Diffusion-toolbox found!?\n\nNB You should *not* ' ...
                     'add the ''old'' sub-directory from the diffusion_toolbox ' ...
                     'directory to your Matlab-path...']))
end
if strcmpi(spm('ver'), 'SPM12')
	disp('=> Adding the old normalization toolbox to your path')
    addpath(fullfile(spm('Dir'),'toolbox','OldNorm'))
end

% Create/update handles structure
handles.output = hObject;	% Choose default command line output for dd_basicproc
handles.DICOM  = [];		% Create userdata variable
handles.Nifti  = [];		% Create userdata variable
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function dd_basicproc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dd_basicproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject, 'Name', 'dd_basicproc 4.0')


% --- Outputs from this function are returned to the command line.
function varargout = dd_basicproc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% ------------------- Data Panel -------------------

% --- Executes on selection change in DataMenu.
function DataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to DataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update data set
update_data_panel(handles, Inf)


% --- Executes during object creation, after setting all properties.
function DataMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddSubjButton.
function AddSubjButton_Callback(hObject, eventdata, handles, PrevSel, PrevT1, Subjects, Series)
% hObject    handle to AddSubjButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adds a row or a column or edits a cell

if nargin<=6
	Series = 1:numel(get(handles.MeasurementMenu, 'String'));
end
if nargin<=5
	Subjects = numel(get(handles.SubjMenu, 'String')) + 1;	% Append a subject row
end
if nargin<=3
	PrevT1  = [];
	PrevSel = [];
end

% Get datatype
DataTypes = get(handles.DataMenu, 'String');
Selected  = get(handles.DataMenu, 'Value');
WDir      = get(handles.DirText, 'String');
if isempty(WDir)
    WDir = pwd;
end

% Get the files & BET-settings from the user
for SubjNr = Subjects
	for SeriesNr = Series
		switch DataTypes{Selected}
			case 'DICOM'
				[FList Sts] = spm_select(Inf, 'any', 'Select DTI data or Job-file (-> edit filter)', ...
										 PrevSel, WDir, '(?i)\.ima$|(?i)\.dcm$');
			case 'Nifti'
				[FList Sts] = spm_select(Inf, 'any', 'Select DTI data or Job-file (-> edit filter)', ...
										 PrevSel, WDir, '\.img$|\.nii$');
		end
		if Sts==false, return, end
		[FPath, ~, Ext] = fileparts(FList(1,:));
		if strcmp(Ext, '.mat')
			for n = 1:size(FList,1)
				[handles Subjects] = appendjob(handles, FList(n,:));
			end
			break
		else
			handles.(DataTypes{Selected})(SubjNr, SeriesNr).Path  = FPath;
			handles.(DataTypes{Selected})(SubjNr, SeriesNr).Files = makefcell(FList);
			if get(handles.T1Box, 'Value')
				switch DataTypes{Selected}
					case 'DICOM'
						T1Text = spm_select(Inf, 'any', 'Select T1 reference scan (DICOM or Nifti)', ...
											PrevT1, FPath,'.nii$|.img$|(?i)\.ima$|(?i)\.dcm$');
					case 'Nifti'
						T1Text = spm_select(1, 'image', 'Select T1 reference scan', PrevT1, FPath);
				end
			else
				T1Text = '';
			end
			handles.((DataTypes{Selected}))(SubjNr, SeriesNr).T1Text = char(strtok(cellstr(T1Text), ','));  % Throw away frame-info
		end
		if get(handles.BETBox, 'Value')
			MenuStr = get(handles.BETMenu, 'String');
			MenuSel = questdlg('What image would you like to use to make a mask?','BET mask', MenuStr{:}, MenuStr{get(handles.BETMenu, 'Value')});
			handles.BETMenu_.Val(SubjNr,SeriesNr) = find(strcmp(MenuSel, MenuStr));
			BETOpts = inputdlg('BET arguments', 'BET mask', 1, {get(handles.BETOpts, 'String')});
			if isempty(BETOpts)
				handles.BETOpts_.Str{SubjNr,SeriesNr} = get(handles.BETOpts, 'String');
			else
				handles.BETOpts_.Str{SubjNr,SeriesNr} = char(BETOpts);
			end
		else
			handles.BETMenu_.Val(SubjNr,SeriesNr) = get(handles.BETMenu, 'Value');
			handles.BETOpts_.Str{SubjNr,SeriesNr} = get(handles.BETOpts, 'String');
		end
	end
end

% Store the data
guidata(hObject, handles)

% Update the gui
update_data_panel(handles, Subjects(end), Series(1))
update_preprocessing_panel(handles)
update_longitudinal_panel(handles)


% --- Executes on button press in DelSubjButton.
function DelSubjButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelSubjButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(get(handles.SubjMenu, 'String'))
	return
end

% Get datatype and subject number
DataTypes	= get(handles.DataMenu, 'String');
Selected	= get(handles.DataMenu, 'Value');
NotSelected = 3 - Selected;		% Assumes Selected is either 1 or 2
SubjNr		= get(handles.SubjMenu, 'Value');

% Delete the current subject
handles.(DataTypes{Selected})(SubjNr,:) = [];
try
	handles.(DataTypes{NotSelected})(SubjNr,:) = [];
	if numel(handles.DICOM) ~= numel(handles.Nifti)
		warning('DIDI:GUI:Conflict', 'DICOM and Nifti subject lists are incompatible')
	end
end
handles.BETMenu_.Val(SubjNr,:) = [];
handles.BETOpts_.Str(SubjNr,:) = [];

% Store the data
guidata(hObject, handles)

% Update the gui
update_data_panel(handles, Inf)
update_preprocessing_panel(handles)


% --- Executes on selection change in SubjMenu.
function SubjMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SubjMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update data set
update_data_panel(handles)
update_preprocessing_panel(handles)


% --- Executes during object creation, after setting all properties.
function SubjMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubjMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DirText_Callback(hObject, eventdata, handles)
% hObject    handle to DirText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get datatype
DataTypes = get(handles.DataMenu, 'String');
Selected  = get(handles.DataMenu, 'Value');
SubjNr    = get(handles.SubjMenu, 'Value');
SeriesNr  = get(handles.MeasurementMenu, 'Value');

% Get the data from the user
handles.(DataTypes{Selected})(SubjNr, SeriesNr).Path = get(hObject, 'String');

% Store the data
guidata(hObject, handles)


function FileText_Callback(hObject, eventdata, handles)
% hObject    handle to FileText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the existing file data
WDir    = get(handles.DirText, 'String');
Files   = get(handles.FileText, 'String');
if isempty(Files)
    return
end
PrevSel = prepend(Files, [WDir filesep]);

% And the T1 data
DataTypes = get(handles.DataMenu, 'String');
Selected  = get(handles.DataMenu, 'Value');
SubjNr    = get(handles.SubjMenu, 'Value');
SeriesNr  = get(handles.MeasurementMenu, 'Value');
switch DataTypes{Selected}
    case 'DICOM'
        PrevT1 = cellstr(handles.DICOM(SubjNr, SeriesNr).T1Text);
    case 'Nifti'
        PrevT1 = {handles.Nifti(SubjNr, SeriesNr).T1Text};
    otherwise
        error('Data-type option not recognized')
end

% And pass it on
AddSubjButton_Callback(hObject, eventdata, handles, PrevSel, PrevT1, ...
	get(handles.SubjMenu, 'Value'), get(handles.MeasurementMenu, 'Value'))


% --- Executes on button press in T1Box.
function T1Box_Callback(hObject, eventdata, handles)
% hObject    handle to T1Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
	OK = check_t1(hObject, handles);
	if isempty(OK)
		set(hObject, 'Value', 0)
	end
end

update_longitudinal_panel(handles)
update_preprocessing_panel(handles)


function T1Text_Callback(hObject, eventdata, handles)
% hObject    handle to T1Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get datatype
DataTypes = get(handles.DataMenu, 'String');
Selected  = get(handles.DataMenu, 'Value');
SubjNr    = get(handles.SubjMenu, 'Value');
SeriesNr  = get(handles.MeasurementMenu, 'Value');

% Get the data from the user
switch DataTypes{Selected}
    case 'DICOM'
        if size(handles.DICOM(SubjNr, SeriesNr).T1Text) > 1
            handles.DICOM(SubjNr, SeriesNr).T1Text = spm_select(Inf, 'any', ...
                'Select T1 reference scan (DICOM or Nifti)', {handles.DICOM(SubjNr, SeriesNr).T1Text}, ...
                fileparts(get(hObject, 'String')), '\.nii$|\.img$|(?i)\.ima$|(?i)\.dcm$');
            update_data_panel(handles)
        else
            handles.DICOM(SubjNr, SeriesNr).T1Text = get(hObject, 'String');
        end
    case 'Nifti'
        handles.Nifti(SubjNr, SeriesNr).T1Text = get(hObject, 'String');
    otherwise
        error('Data-type option not recognized')
end

% Store the data
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function T1Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T1Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% ------------------- Longitudinal Panel -------------------

% --- Executes on button press in AddMeasureButton.
function AddMeasureButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddMeasureButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SubjNr	 = get(handles.SubjMenu, 'Value');
SeriesNr = numel(get(handles.MeasurementMenu,'String')) + 1;
if isempty(get(handles.SubjMenu, 'String'))
	SeriesNr = 1:2;			% Catch state where Subject 1, Series 1 is empty (e.g. vanilla start)
end

AddSubjButton_Callback(hObject, eventdata, handles, [], [], SubjNr, SeriesNr)


% --- Executes on button press in DelMeasureButton.
function DelMeasureButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelMeasureButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if numel(get(handles.MeasurementMenu, 'String'))== 1
	warndlg('You cannot have less than one measurement', 'Remove measurement')
	return
end

% Get datatype and series number
DataTypes = get(handles.DataMenu, 'String');
Selected  = get(handles.DataMenu, 'Value');
NotSelected = 3 - Selected;		% Assumes Selected is either 1 or 2
SeriesNr  = get(handles.MeasurementMenu, 'Value');

% Delete the current series
handles.(DataTypes{Selected})(:,SeriesNr) = [];
try
	handles.(DataTypes{NotSelected})(:,SeriesNr) = [];
	if numel(handles.DICOM) ~= numel(handles.Nifti)
		warning('DIDI:GUI:Conflict', 'DICOM and Nifti subject lists are incompatible')
	end
end
handles.BETMenu_.Val(:,SeriesNr) = [];
handles.BETOpts_.Str(:,SeriesNr) = [];

% Store the data
guidata(hObject, handles)

% Update the gui
update_data_panel(handles, [], Inf)
update_preprocessing_panel(handles)
update_longitudinal_panel(handles)


% --- Executes on selection change in MeasurementMenu.
function MeasurementMenu_Callback(hObject, eventdata, handles)
% hObject    handle to MeasurementMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MeasurementMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MeasurementMenu

update_data_panel(handles)
update_preprocessing_panel(handles)


% --- Executes during object creation, after setting all properties.
function MeasurementMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasurementMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in T1RealignBox.
function T1RealignBox_Callback(hObject, eventdata, handles)
% hObject    handle to T1RealignBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T1RealignBox


% --- Executes on button press in SinglePEUnwarpBox.
function SinglePEUnwarpBox_Callback(hObject, eventdata, handles)
% hObject    handle to SinglePEUnwarpBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SinglePEUnwarpBox


% --- Executes on selection change in HyperalignMenu.
function HyperalignMenu_Callback(hObject, eventdata, handles)
% hObject    handle to HyperalignMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HyperalignMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HyperalignMenu


% --- Executes during object creation, after setting all properties.
function HyperalignMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HyperalignMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% ------------------- Processing Panel -------------------


% --- Executes on selection change in DenoisingMenu.
function DenoisingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to DenoisingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_preprocessing_panel(handles)


% --- Executes during object creation, after setting all properties.
function DenoisingMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DenoisingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RicianBox.
function RicianBox_Callback(hObject, eventdata, handles)
% hObject    handle to RicianBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RicianBox


% --- Executes on selection change in BETMenu.
function BETMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BETMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BETMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BETMenu

% Get the data from the user
if get(handles.BETBox, 'Value')
	SubjNr   = get(handles.SubjMenu, 'Value');
	SeriesNr = get(handles.MeasurementMenu, 'Value');
else
	SubjNr   = 1:numel(get(handles.SubjMenu, 'String'));
	SeriesNr = 1:numel(get(handles.MeasurementMenu, 'String'));
end	
for n = SubjNr
	for m = SeriesNr
		handles.BETMenu_.Val(n,m) = get(hObject, 'Value');
	end
end

% Store the data
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function BETMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BETMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Get the defaults data from the GUI
handles.BETMenu_.Val = get(hObject, 'Value');

% Store the data
guidata(hObject, handles)


function BETOpts_Callback(hObject, eventdata, handles)
% hObject    handle to BETOpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BETOpts as text
%        str2double(get(hObject,'String')) returns contents of BETOpts as a double

% Get the data from the user
if get(handles.BETBox, 'Value')
	SubjNr   = get(handles.SubjMenu, 'Value');
	SeriesNr = get(handles.MeasurementMenu, 'Value');
else
	SubjNr   = 1:numel(get(handles.SubjMenu, 'String'));
	SeriesNr = 1:numel(get(handles.MeasurementMenu, 'String'));
end	
for n = SubjNr
	for m = SeriesNr
		handles.BETOpts_.Str{n,m} = get(hObject, 'String');
	end
end

% Store the data
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function BETOpts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BETOpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Get the defaults data from the GUI
handles.BETOpts_.Str{1,1} = get(hObject, 'String');

% Store the data
guidata(hObject, handles)


% --- Executes on button press in BETBox.
function BETBox_Callback(hObject, eventdata, handles)
% hObject    handle to BETBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
	return
else
	SubjNr   = 1:numel(get(handles.SubjMenu, 'String'));
	SeriesNr = 1:numel(get(handles.MeasurementMenu, 'String'));
end	
for n = SubjNr
	for m = SeriesNr
		handles.BETMenu_.Val(n,m) = get(handles.BETMenu, 'Value');
		handles.BETOpts_.Str{n,m} = get(handles.BETOpts, 'String');
	end
end

% Store the data
guidata(hObject, handles)

update_preprocessing_panel(handles)


% --- Executes on selection change in ArtDetMenu.
function ArtDetMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ArtDetMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_preprocessing_panel(handles)
update_tensor_panel(handles)


% --- Executes during object creation, after setting all properties.
function ArtDetMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ArtDetMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PATCHText_Callback(hObject, eventdata, handles)
% hObject    handle to PATCHText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PATCHText as text
%        str2double(get(hObject,'String')) returns contents of PATCHText as a double


% --- Executes during object creation, after setting all properties.
function PATCHText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PATCHText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RealignMenu.
function RealignMenu_Callback(hObject, eventdata, handles)
% hObject    handle to RealignMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function RealignMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RealignMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PEUnwarpBox.
function PEUnwarpBox_Callback(hObject, eventdata, handles)
% hObject    handle to PEUnwarpBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_longitudinal_panel(handles)
update_preprocessing_panel(handles)


function PEUnwarpText_Callback(hObject, eventdata, handles)
% hObject    handle to PEUnwarpText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PEUnwarpText as text
%        str2double(get(hObject,'String')) returns contents of PEUnwarpText as a double

if numel(str2num(get(hObject,'String')))~=3
	warndlg('Order must be a three element numerical array')
	set(hObject,'String','9 9 9')
end


% --- Executes during object creation, after setting all properties.
function PEUnwarpText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PEUnwarpText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CaminoBox.
function CaminoBox_Callback(hObject, eventdata, handles)
% hObject    handle to CaminoBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CaminoBox


% --- Executes on button press in FSLBox.
function FSLBox_Callback(hObject, eventdata, handles)
% hObject    handle to FSLBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% ------------------- Tensor Estimation Panel -------------------

% --- Executes on selection change in EstMenu.
function EstMenu_Callback(hObject, eventdata, handles)
% hObject    handle to EstMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EstMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EstMenu
% If 'none' is chosen, check if PATCH estimation is used

update_tensor_panel(handles)
update_preprocessing_panel(handles)


% --- Executes during object creation, after setting all properties.
function EstMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EstMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MaxBValText_Callback(hObject, eventdata, handles)
% hObject    handle to MaxBValHeader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxBValHeader as text
%        str2double(get(hObject,'String')) returns contents of MaxBValHeader as a double


% --- Executes during object creation, after setting all properties.
function MaxBValText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxBValHeader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MDBox.
function MDBox_Callback(hObject, eventdata, handles)
% hObject    handle to MDBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in NormDBox.
function NormDBox_Callback(hObject, eventdata, handles)
% hObject    handle to NormDBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NormDBox


% --- Executes on button press in NormABox.
function NormABox_Callback(hObject, eventdata, handles)
% hObject    handle to NormABox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NormABox


% --- Executes on button press in FABox.
function FABox_Callback(hObject, eventdata, handles)
% hObject    handle to FABox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RABox.
function RABox_Callback(hObject, eventdata, handles)
% hObject    handle to RABox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ModeBox.
function ModeBox_Callback(hObject, eventdata, handles)
% hObject    handle to ModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ModeBox


% --- Executes on button press in EigBox.
function EigBox_Callback(hObject, eventdata, handles)
% hObject    handle to EigBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(hObject,'Value')
	set(handles.RDBox, 'Value', 0)
end


% --- Executes on button press in RDBox.
function RDBox_Callback(hObject, eventdata, handles)
% hObject    handle to RDBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
	set(handles.EigBox, 'Value', 1)
end


% --- Executes on button press in NormBox.
function NormBox_Callback(hObject, eventdata, handles)
% hObject    handle to NormBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NormBox


%% ------------------- Job Panel -------------------

% --- Executes on button press in CleanUpBox.
function CleanUpBox_Callback(hObject, eventdata, handles)
% hObject    handle to CleanUpBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CleanUpBox


% --- Executes on button press in ParallelBox.
function ParallelBox_Callback(hObject, eventdata, handles)
% hObject    handle to ParallelBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ParallelBox
if ~exist('qsubcellfun', 'file')
	addpath(fullfile(filesep,'home','common','matlab','fieldtrip','qsub'), '-end')
	if ~exist('qsubcellfun', 'file')
		warndlg('No distributed computing functions found. Make sure you added these to your matlabpath.')
		set(hObject, 'Value', 0)
		return
	end
end


% --- Executes on button press in HelpButton.
function HelpButton_Callback(hObject, eventdata, handles)
% hObject    handle to HelpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

doc(mfilename)


% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles, Job)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<4 || isempty(Job)
	[FName Path] = uigetfile({'*.mat'}, 'Load previous job settings', ['dd_' datestr(now,'yyyymmmdd') '_Job.mat']);
	if ~FName, return, end
	load(fullfile(Path, FName))
end
if ~strcmp(get(handles.dd_basicproc_gui, 'Name'), Job.Ver)
    Ans = questdlg(['Another version (' Job.Ver ') was used to create this Job. Use the settings anyway (errors may occur)?']);
    if ~strcmp(Ans, 'Yes'), return, end
end

% -- Data panel --
handles.DICOM = Job.DICOM;
handles.Nifti = Job.Nifti;
tryset(handles,	Job, 'DataMenu')
tryset(handles, Job, 'T1Box')
% -- Longitudinal panel --
tryset(handles, Job, 'T1RealignBox')
tryset(handles, Job, 'SinglePEUnwarpBox');
tryset(handles, Job, 'HyperalignMenu');
% -- Pre-processing panel --
tryset(handles, Job, 'DenoisingMenu')
tryset(handles, Job, 'RicianBox')
if numel(get(handles.BETMenu,'String'))==numel(Job.BETMenu.Str) && all(strcmpi(get(handles.BETMenu,'String'), Job.BETMenu.Str))
	handles.BETMenu_.Val = Job.BETMenu.Val;
	handles.BETOpts_.Str = cellstr(Job.BETOpts.Str);
else
	warning('Could not set job-items: BETMenu & BETOpts')
	handles.BETMenu_.Val = get(handles.BETMenu,'Value');
	handles.BETOpts_.Str = {get(handles.BETOpts,'String')};
end
NrSubjSeries = max([size(Job.DICOM);size(Job.Nifti)]);
if numel(handles.BETMenu_.Val)==1 && prod(NrSubjSeries)>1
	for SubjNr = 1:NrSubjSeries(1)
		for SeriesNr = 1:NrSubjSeries(2)
			handles.BETMenu_.Val(SubjNr,SeriesNr) = get(handles.BETMenu,'Value');
			handles.BETOpts_.Str{SubjNr,SeriesNr} = get(handles.BETOpts,'String');
		end
	end
end
tryset(handles, Job, 'BETBox')
tryset(handles, Job, 'ArtDetMenu')
tryset(handles, Job, 'PATCHText')
tryset(handles, Job, 'RealignMenu')
tryset(handles, Job, 'PEUnwarpBox')
tryset(handles, Job, 'PEUnwarpText')
tryset(handles, Job, 'CaminoBox')
tryset(handles, Job, 'FSLBox')
% -- Tensor panel --
tryset(handles, Job, 'EstMenu')
tryset(handles, Job, 'MaxBValText')
tryset(handles, Job, 'MDBox')
tryset(handles, Job, 'NormDBox')
tryset(handles, Job, 'NormABox')
tryset(handles, Job, 'FABox')
tryset(handles, Job, 'RABox')
tryset(handles, Job, 'ModeBox')
tryset(handles, Job, 'EigBox')
tryset(handles, Job, 'RDBox')
tryset(handles, Job, 'NormBox')
% -- Job panel --
tryset(handles, Job, 'CleanUpBox')
tryset(handles, Job, 'ParallelBox')

guidata(hObject, handles)

update_data_panel(handles, Inf, Inf)
update_longitudinal_panel(handles)
update_preprocessing_panel(handles)
update_tensor_panel(handles)


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do some input checking
onCleanup(@()set(hObject, 'Enable', 'on', 'ForegroundColor', [0 0 0]));
if get(handles.T1Box, 'Value')
    OK = check_t1(hObject, handles);
    if isempty(OK) || ~OK
        return
    end
end
MenuStr = get(handles.DataMenu, 'String');
if strcmpi(MenuStr{get(handles.DataMenu, 'Value')}, 'Nifti')
    
	if any(strncmp('artc_', cat(1,handles.Nifti.Files), 5))
		switch questdlg('Would you like to use the artefact-corrected (artc_*) files?')
			case 'No'
				for n = 1:numel(handles.Nifti)
					handles.Nifti(n).Files = strrep(handles.Nifti(n).Files, 'artc_','');
				end
				guidata(hObject, handles)
				update_data_panel(handles)
			case {'Cancel' ''}
				return
		end
	end
    
	set(hObject, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
	HW = waitbar(0, 'Checking volume orientations...');
	for n = 1:numel(handles.Nifti)
		WDir = [handles.Nifti(n).Path filesep];
		if numel(handles.Nifti(n).Files) == 1 && ~exist(fullfile(WDir, char(handles.Nifti(n).Files)), 'file')	% Catch cases where the 4D file has been split but the job was not updated (e.g. due to a crash)
			fprintf('DW Image %s not found, using splitted 3D images instead\n', handles.Nifti(n).Files{:})
			handles.Nifti(n).Files = cellstr(spm_select('List', WDir, [char(spm_file(handles.Nifti(n).Files,'basename')) '_\d\d\d\d\d\.nii$']));
			guidata(hObject, handles)
			update_data_panel(handles)
		end
		FList = prepend(char(handles.Nifti(n).Files), WDir);
		try
			spm_check_orientations(spm_vol(FList));
		catch
			fprintf(['One or more DW images do not have the same orientation\n' ...
					 'Resetting orientations in: %s\n'], WDir)
			try
				for m = 1:size(FList,1)
					D = dti_get_dtidata(FList(m,:));
					if isfield(D,'mat')
						spm_get_space(FList(m,:), D.mat);
					elseif isfield(D,'M')		% Convert my earlier M-convention
						spm_get_space(FList(m,:), D.M);
						D.mat = D.M;
						dti_get_dtidata(FList(m,:), rmfield(D,'M'));
					end
				end
				spm_check_orientations(spm_vol(FList));
			catch
				warndlg('Orientations are messed up. Use DICOM images instead')
				myclose(HW)
				set(hObject, 'Enable', 'on', 'ForegroundColor', [0 0 0])
 				return
			end
		end
		mywaitbar(n/numel(handles.Nifti), HW)
	end
	myclose(HW)
	set(hObject, 'Enable', 'on', 'ForegroundColor', [0 0 0])
    
end

% Create the job structure
% -- Data panel --
Job.Ver						= get(handles.dd_basicproc_gui, 'Name');
Job.DICOM					= handles.DICOM;
Job.Nifti					= handles.Nifti;
Job.DataMenu.Val			= get(handles.DataMenu, 'Value');
Job.DataMenu.Str			= get(handles.DataMenu, 'String');
Job.T1Box.Val				= get(handles.T1Box, 'Value');
% -- Longitudinal panel --
Job.T1RealignBox.Val		= get(handles.T1RealignBox, 'Value');
Job.SinglePEUnwarpBox.Val	= get(handles.SinglePEUnwarpBox, 'Value');
Job.HyperalignMenu.Val		= get(handles.HyperalignMenu, 'Value');
Job.HyperalignMenu.Str		= get(handles.HyperalignMenu, 'String');
% -- Pre-processing panel --
Job.DenoisingMenu.Val		= get(handles.DenoisingMenu, 'Value');
Job.DenoisingMenu.Str		= get(handles.DenoisingMenu, 'String');
Job.RicianBox.Val			= get(handles.RicianBox, 'Value');
Job.BETMenu.Val				= handles.BETMenu_.Val;
Job.BETMenu.Str				= get(handles.BETMenu, 'String');
Job.BETOpts.Str				= handles.BETOpts_.Str;
Job.BETBox.Val				= get(handles.BETBox, 'Value');
Job.ArtDetMenu.Val			= get(handles.ArtDetMenu, 'Value');
Job.ArtDetMenu.Str			= get(handles.ArtDetMenu, 'String');
Job.PATCHText.Str			= get(handles.PATCHText, 'String');
Job.RealignMenu.Val			= get(handles.RealignMenu, 'Value');
Job.RealignMenu.Str			= get(handles.RealignMenu, 'String');
Job.PEUnwarpBox.Val			= get(handles.PEUnwarpBox, 'Value');
Job.PEUnwarpText.Str		= get(handles.PEUnwarpText, 'String');
Job.CaminoBox.Val			= get(handles.CaminoBox, 'Value');
Job.FSLBox.Val				= get(handles.FSLBox, 'Value');
% -- Tensor panel --
Job.EstMenu.Val				= get(handles.EstMenu, 'Value');
Job.EstMenu.Str				= get(handles.EstMenu, 'String');
Job.MaxBValText.Str			= get(handles.MaxBValText, 'String');
Job.MDBox.Val				= get(handles.MDBox, 'Value');
Job.NormDBox.Val			= get(handles.NormDBox, 'Value');
Job.NormABox.Val			= get(handles.NormABox, 'Value');
Job.FABox.Val				= get(handles.FABox, 'Value');
Job.RABox.Val				= get(handles.RABox, 'Value');
Job.ModeBox.Val				= get(handles.ModeBox, 'Value');
Job.EigBox.Val				= get(handles.EigBox, 'Value');
Job.RDBox.Val				= get(handles.RDBox, 'Value');
Job.NormBox.Val				= get(handles.NormBox, 'Value');
% -- Job panel --
Job.CleanUpBox.Val			= get(handles.CleanUpBox, 'Value');
Job.ParallelBox.Val			= get(handles.ParallelBox, 'Value');

% Get the name of the file to store the Job-settings in (see jobman)
[FName Path] = uiputfile({'*.mat'}, 'Save job settings', ['dd_' datestr(now,'yyyymmmdd') '_Job']);
if ~FName
    return
end

% Check if we have a distributed computing infrastructure that can do the job
if Job.ParallelBox.Val && ~exist('qsubcellfun', 'file')
	addpath(fullfile(filesep,'home','common','matlab','fieldtrip','qsub'), '-end')
	if ~exist('qsubcellfun', 'file')
		warndlg('No parallel computing functions found. Make sure you added these to your matlabpath.')
		return
	end
end

% Run the job!
set(hObject, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
[Failed Job] = jobman(Job, fullfile(Path, FName));

% Update the GUI
if ishandle(hObject)
	handles.Nifti = Job.Nifti;
	guidata(hObject, handles)
	update_data_panel(handles, Inf)
	set(hObject, 'Enable', 'on', 'ForegroundColor', [0 0 0])
end


%% ------------------- DTI Job -------------------

function [Failed Job] = jobman(Job, JobName)

if nargin<2 || isempty(JobName)
	JobName = fullfile(pwd, ['dd_' datestr(now,'yyyymmmdd') '_Job.mat']);
end

% Store the job info on disk and define the LogNames
[Path FName] = fileparts(JobName);						% Discard the file-extension
JobName		 = fullfile(Path, [FName '.mat']);			% And make sure the extension is .mat
LogName0	 = fullfile(Path, [FName '_LogS0000.ps']);
NrSubj		 = size(Job.(Job.DataMenu.Str{Job.DataMenu.Val}), 1);
for n = 1:NrSubj
	LogNames{n}	= sprintf('%s%04d.ps', LogName0(1:end-7), n);
	% Clean-up previous log-files
	if exist(LogNames{n}, 'file')
		delete(LogNames{n})
	end
	if exist([LogNames{n}(1:end-2) 'tsv'], 'file')
		delete([LogNames{n}(1:end-2) 'tsv'])
	end
end
disp(['Saving Job structure in: ' JobName])
save(JobName, 'Job')

% Run the job!
if ishandle(gcbf)
	print(gcbf, [LogName0(1:end-2) 'pdf'], '-dpdf', '-loose')
end
DiaryFile = fullfile(Path, [FName '_diary.txt']);
if exist(DiaryFile, 'file')
	delete(DiaryFile)
end
diary(DiaryFile)
fprintf('\n\n************ BEGIN DD_BASICPROC JOB (%s) ************\n\n', datestr(now))
t0 = clock;
HG = spm_figure('GetWin', 'Graphics'); set(HG, 'Visible', 'Off')	% Create (invisible) graphics windows
HI = spm_figure('GetWin', 'Interactive'); set(HI, 'Visible', 'Off')
% onCleanup(@()close(HG,HI));
try
	spm_get_defaults;				% Does not work with spm5
catch
	global defaults;
	if isempty(defaults)
		spm_defaults;
	end
end
if Job.ParallelBox.Val
	spm_get_defaults('cmdline',true)
end
spm_jobman('initcfg');						% Setup for batch system
[RunJob TDImgs] = run_job(Job, LogNames);	% RunJob can be used to replace the vannila nii-files in the saved Job with filtered or corrected versions (currently not used)
if strcmp(Job.DataMenu.Str{Job.DataMenu.Val}, 'DICOM')
	% We only did DICOM conversion and returned here to save the info to disk
	Job.Nifti        = RunJob.Nifti;
	Job.DataMenu.Val = find(strcmp('Nifti', Job.DataMenu.Str));  % Use the Nifti-files next time
	save(JobName, 'Job')
	% Now restart the job
	[RunJob TDImgs] = run_job(Job, LogNames);
end

% Save the names of the (warped) tensor-derivative images for subsequent analyses
Failed = cellfun(@isempty,{Job.Nifti.Files});
if ~strcmp(Job.EstMenu.Str{Job.EstMenu.Val}, 'none')
	TDVal	= logical([Job.MDBox.Val Job.NormDBox.Val Job.NormABox.Val Job.FABox.Val Job.RABox.Val Job.ModeBox.Val Job.RDBox.Val Job.RDBox.Val]);	% NB: Match TDVal with DD_BASICPROC_ESTDTENSOR
	PostFix = {'_ad' '_nd' '_na' '_fa' '_va' '_mo' '_eval1' '_eval23'};		% NB: Match order with DD_BASICPROC_ESTDTENSOR
	PostFix = PostFix(TDVal);
	for n = 1:numel(PostFix)
		fprintf('\n==> Saving the names of the %s images for subsequent analyses', PostFix{n})
		FID = fopen([fullfile(Path, FName) PostFix{n} '.txt'], 'w');
		for m = 1:length(TDImgs)
			if isempty(TDImgs{m})
				Failed(m) = true;
			else
				fprintf(FID, '%s\n', TDImgs{m}(n,:));
			end
		end
		fclose(FID);
	end
end
Failed = any(reshape(Failed,[],size(Job.Nifti,2)), 2);

% Convert all postscript logfiles to pdf and delete the ps-file if we have a pdf
fprintf('\n\n==> Converting all postscript logfiles to pdf\n')
for n = 1:numel(LogNames)
	if exist(LogNames{n},'file') && ~system(['ps2pdf "' LogNames{n} '" "' LogNames{n}(1:end-2) 'pdf"'])
		delete(LogNames{n})
		fprintf('.')
	else
		fprintf('!')
	end
end

% Update the Job with the new Nifti- and Output-info and save it to disk
if any(Failed)					% Package and save failed jobs for reprocessing
	warning('DIDI:FailedSubjects', ['\nDistributed processing failed to complete in %d dataset(s):\n%s. You can use ' ...
									'%s_failed.mat jobfile for a retry...'], sum(Failed), ...
									sprintf('%s\n',Job.Nifti(Failed,:).Path), JobName(1:end-4))
	FailedJob = Job;
	FinishJob = Job;
	if ~isempty(Job.DICOM)
		FailedJob.DICOM = Job.DICOM(Failed,:);
		FinishJob.DICOM = Job.DICOM(~Failed,:);
	end
	FailedJob.Nifti = Job.Nifti(Failed,:);
	FinishJob.Nifti = Job.Nifti(~Failed,:);
	if ~Job.CleanUpBox.Val
		FinishJob.Output = RunJob.Output(~Failed,:);
	end
	Tmp = Job;
	Job = FailedJob; save([JobName(1:end-4) '_failed'], 'Job')
	Job = FinishJob; save([JobName(1:end-4) '_successful'], 'Job')
	Job = Tmp;
elseif ~Job.CleanUpBox.Val			% Everything went allright and we have intermediate files
	Job.Output = RunJob.Output;
	save(JobName, 'Job')
end

% We're done :-)
ETime = etime(clock, t0);
fprintf('\n************ BEGAN DD_BASICPROC JOB (%s) ************', datestr(t0))
fprintf('\n************ END DD_BASICPROC JOB (%s)   ************', datestr(now))
fprintf('\n************ ELAPSED TIME: %02g hour and %2.2f minutes\t   ************\n', ...
		floor(ETime/3600), rem(ETime,3600)/60)
diary('Off')


function [Job TDImgs] = run_job(Job, LogNames)

% All imaging sessions of a single subject are processed as a set. The pipeline
% has the following stages:
%
% 1) DICOM => Nifti
%    - file conversion
%    - return to caller to first save the Job-info (in case stage 2 crashes)
% 2) Nifti => rNifti
%	 - realign the T1 series
%    - SNR & b0-mask
%	 - denoising
%    - realignment
%    - artefact detection
% 3) rNifti => tensor estimation (= fully automatic)
%    - realignment
%    - normalization
%    - compute Dxx, Dxy, ...
%    - compute tensor derivatives
%    - compute eigenvectors & values
%	 - hyper-align and normalize the tensor derivative images

% Set (sq DT-error) defaults and open graphical output windows
Constr.sptl.apply  = true;		% Use spatial constraints
Constr.sptl.basset = 'cos';		% Basis set ('cos'/'pol')
%Constr.sptl.order  = 4;		% 4th order spatial basis function
Constr.sptl.order  = 2;			% 2nd order spatial basis function
Constr.grad.apply  = true;		% Use gradient constraints
%Constr.grad.shear  = [1 2 3];	% Linear effects + cross-terms
%Constr.grad.scale  = [1 2 3];	% Linear effects + cross-terms
%Constr.grad.trans  = [1 2 3];	% Linear effects + cross-terms
Constr.grad.shear  = 1;			% Linear effects only
Constr.grad.scale  = 1;			% Linear effects only
Constr.grad.trans  = 1;			% Linear effects only
Tuning.fwhm        = 10;		% Smoothing kernel (default 4/15)
Tuning.imask.apply = true;		% Use "importance masking"
Tuning.imask.type  = 1;			% A clever but slower mask (2 = maxgradient and 3 = maxintensity)
Tuning.imask.prcnt = 15;		% < 20 discards FOV crap (defaults 10/50)
Tuning.difm        = true;		% Base R on diffusion tensor model
Tuning.orth        = [];		% Effects not be included in R
HG = spm_figure('FindWin', 'Graphics');
HI = spm_figure('FindWin', 'Interactive');
spm_figure('Clear', HI);
spm('FigName', 'Processing data', HI);
HW	= [];
HWS = [];
if ishandle(HI)					% Create some more interactive windows (waitbars)
	Units = get(HI, 'Units');
	myset(HG, 'Visible', 'Off')
	myset(HI, 'Units', 'Pixels', 'Visible', 'Off')
	WaitPos = get(HI, 'Position');	% Get the position in pixels
	myset(HI, 'Units', Units)
	WaitPos	= [WaitPos(1:3) 75] + [0 WaitPos(4)+65 0 0];
	if (numel(Job.DICOM) > 1 || numel(Job.Nifti) > 1) &&  ~Job.ParallelBox.Val
		HW = waitbar(0, ['Converting DICOM files (' num2str(numel(Job.DICOM)) ' data sets) ...'], ...
			'Name','dd_basicproc progress-bar', 'Units', 'Pixels', 'Position', WaitPos);
		WaitPos = WaitPos + [0 WaitPos(4)+30 0 0];
	end
	HWS = waitbar(0, 'Processing...', 'Name','Sub-progress bar', 'Units', 'Pixels', ...
		'Position', WaitPos, 'Visible', 'Off');
end


%% --- Processing: Stage 1 (DICOM Job) ---

if strcmp(Job.DataMenu.Str{Job.DataMenu.Val}, 'DICOM')

    % Delete unselected job info
    Job.Nifti = [];

	% DICOM => Nifti file conversion
	fprintf('\n==> Converting DICOM to Nifti (%s)\n', datestr(now))
	if Job.ParallelBox.Val
		
		% Determine the maximum number of DW directions (images)
		NDir   = max(cellfun(@numel,{Job.DICOM(:).Files}));
		[TimReq MemReq] = maxreq(NDir*10, NDir*2*1024^2);
		fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
		[NiList T1Text] = qsubcellfun(@dd_basicproc_convert2nifti, ...
									  stackjob(Job,'DICOM','LinInd'), num2cell(1:numel(Job.DICOM)), ...
									  'TimReq',TimReq, 'MemReq',MemReq, 'jvm','no', ...
									  'StopOnError',false, 'options','-l gres=bandwidth:999');		% Use a linear index for ease (we have no log-output, Job.Nifti will be reshaped later)
		
		% Prepare the Nifti job
		for n = 1:numel(Job.DICOM)
			if ~isempty(NiList{n})
				Job.Nifti(n).Path   = fileparts(NiList{n}(1,:));
				Job.Nifti(n).Files  = makefcell(NiList{n});
				Job.Nifti(n).T1Text = T1Text{n};
			else
				fprintf('Failed: %s\n', Job.DICOM(n).Path)
				Job.Nifti(n).Path   = '';
				Job.Nifti(n).Files  = {};
				Job.Nifti(n).T1Text = '';
			end
		end

	else
		
		myset(HI, 'Visible', 'On')
		for n = 1:numel(Job.DICOM)
			
			spm('FigName', 'DICOM conversion', HI);
			[NiList T1Text] = dd_basicproc_convert2nifti(Job, n);
			
			% Prepare the Nifti job
			Job.Nifti(n).Path   = fileparts(NiList(1,:));
			Job.Nifti(n).Files  = makefcell(NiList);
			Job.Nifti(n).T1Text = T1Text;
			
			mywaitbar(n/numel(Job.DICOM), HW)
			
		end
		myclose(HW, HWS)
		
	end
	Job.Nifti = reshape(Job.Nifti, size(Job.DICOM));	% -> Restore measurements as rows
	TDImgs = '';
    return					% This is to save the new Job-info to disk before continuing

end


%% --- Processing: Stage 2 (Nifti Job) ---

% Create empty logfiles with provenance info
for n = 1:numel(Job.Nifti)
	FIDLog = fopen([LogNames{n}(1:end-2) 'tsv'], 'w');
	fprintf(FIDLog, 'S%d\tPath:\t%s\n', n, Job.Nifti(n).Path);
	fclose(FIDLog);
end

% Clean-up old mask files
mywaitbar(0, HW, 'Cleaning-up old masks')
for n = 1:numel(Job.Nifti)
	delete(fullfile(Job.Nifti(n).Path, 'mean*_*'))
	mywaitbar(n/numel(Job.Nifti), HW)
end

% Convert 4D nifti files into 3D nifti files
for n = 1:numel(Job.Nifti)
	if numel(Job.Nifti(n).Files) == 1
		Hdr = spm_vol(fullfile(Job.Nifti(n).Path, char(Job.Nifti(n).Files)));
		if length(Hdr) > 1
			disp(['Splitting 4D nifti-file: ' Hdr(1).fname])
			Hdrs			   = spm_file_split(Hdr);
			Job.Nifti(n).Files = spm_file({Hdrs.fname}, 'filename');
			delete(Hdr(1).fname)
		end
	end
end

% Determine the maximum number of DW directions (images) and split the file if it is 4D
NDir			  = max(cellfun(@numel,{Job.Nifti.Files}));
[NrSubj	NrSeries] = size(Job.Nifti);
fprintf('Maximum number of DW directions: (%i)\n', NDir)

% Import the bval- bvec info from user provided files if there are no DICOM generated .mat files available
myset(HG, 'Visible', 'On')
myset(HI, 'Visible', 'On')
mywaitbar(0, HW, ['Verifying diffusion information (' num2str(NrSubj*NrSeries) ' subjects)...'])
for n = 1:NrSubj*NrSeries
	dd_basicproc_bvalvec2mat(Job, n);
	mywaitbar(n/(NrSubj*NrSeries), HW)
end

% Realign the T1 series and use the mean T1-image (i.e. replace the original user-entries in the Job)
if Job.T1Box.Val && NrSeries>1 && Job.T1RealignBox.Val
	fprintf('\n==> Realign the T1 series (%s)\n', datestr(now))
	myset(HG, 'Visible', 'On')
	myset(HI, 'Visible', 'On')
	mywaitbar(0, HW, ['Realignment of the ' num2str(NrSeries) ' T1 measurements'])
	for SubjNr = 1:NrSubj
		T1Hdr = spm_realign(char(Job.Nifti(SubjNr,:).T1Text), struct('quality',1, 'sep',2, 'rtm',1));
		spm_reslice(T1Hdr, struct('which',0, 'mean',1))			% Save a mean-T1 volume
		[Job.Nifti(SubjNr,:).T1Text] = deal(prepend(T1Hdr(1).fname, 'mean'));
		mywaitbar(SubjNr/NrSubj, HW)
	end
end

% Estimate the SNR of the raw b0- & DW-images and create brain masks
fprintf('\n==> Estimating the image SNR (%s)\n', datestr(now))
if Job.ParallelBox.Val
	[TimReq MemReq] = maxreq(NDir * 20 * NrSeries, NDir * 80*1024^2);		% Estimated mem 3.6GB, it required 700MB
	fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
	OK = qsubcellfun(@dd_basicproc_getsnr, stackjob(Job,'Nifti','Subj'), num2cell(1:NrSubj), LogNames, ...
					  'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);
else
	myset(HG, 'Visible', 'On')
	myset(HI, 'Visible', 'On')
	mywaitbar(0, HW, ['Noise level estimation (' num2str(NrSubj) ' subjects)...'])
	for SubjNr = 1:NrSubj
		OK{SubjNr} = dd_basicproc_getsnr(Job, SubjNr, LogNames{SubjNr});
		mywaitbar(SubjNr/NrSubj, HW)
	end
end
for SubjNr = 1:NrSubj
	if isempty(OK{SubjNr}) || any(~OK{SubjNr})	% Discard the entire serie
		for SeriesNr = 1:NrSeries
			fprintf('Failed: %s\n', Job.Nifti(SubjNr,SeriesNr).Path)
			Job.Nifti(SubjNr,SeriesNr).Path	  = '';
			Job.Nifti(SubjNr,SeriesNr).Files  = {};
			Job.Nifti(SubjNr,SeriesNr).T1Text = '';
		end
	end
end

if ~strcmp(Job.DenoisingMenu.Str{Job.DenoisingMenu.Val}, 'none')
	
	% Start the image denoising and recreate the masks using the denoised images
	fprintf('\n==> Denoising DW images (%s)\n', datestr(now))
	NiList = makeflist(Job.Nifti);				% cellarray of size NrSubj x NrSeries
	if Job.ParallelBox.Val
		[TimReq MemReq] = maxreq(NDir * 800, NDir * 150*1024^2);
		fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
		DNFiles = qsubcellfun(@dd_basicproc_denoising, NiList, ...
							  repmat(Job.DenoisingMenu.Str(Job.DenoisingMenu.Val), NrSubj, NrSeries), ...
							  repmat({Job.RicianBox.Val}, NrSubj, NrSeries), ...
							  repmat(LogNames, 1, NrSeries), ...
							  repmat(num2cell(1:NrSeries), NrSubj, 1), ...
							  'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);
		for n = 1:numel(Job.Nifti)
			Job.Nifti(n).Files = makefcell(DNFiles{n});		% Update the Job info
		end
		[TimReq MemReq] = maxreq(NDir * 30 * NrSeries, NDir * 30*1024^2);
		fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
		Mask = qsubcellfun(@dd_basicproc_getmask, stackjob(Job,'Nifti','Subj'), num2cell(1:NrSubj), cell(1,NrSubj), ...
						   LogNames, 'TimReq',TimReq, 'MemReq',MemReq);
	else
		mywaitbar(0, HW, ['Image denoising (' num2str(numel(Job.Nifti)) ' data sets)...'])
		for SubjNr = 1:NrSubj
			for SeriesNr = 1:NrSeries
				DNFiles	= dd_basicproc_denoising(NiList{SubjNr,SeriesNr}, Job.DenoisingMenu.Str{Job.DenoisingMenu.Val}, ...
												 Job.RicianBox.Val, LogNames{SubjNr}, SeriesNr);
				Job.Nifti(SubjNr,SeriesNr).Files = makefcell(DNFiles);		% Update the Job info
				mywaitbar(((SubjNr-1)*NrSeries+SeriesNr)/numel(Job.Nifti), HW)
			end
			Mask{SubjNr} = dd_basicproc_getmask(Job, SubjNr, [], LogNames{SubjNr});
		end
	end
	for SubjNr = 1:NrSubj
		for SeriesNr = 1:numel(Mask{SubjNr}) * any(isempty(Mask{SubjNr}))	% Discard the entire serie, otherwise = find(isempty(Mask{SubjNr}))
			fprintf('Failed: %s\n', Job.Nifti(SubjNr,SeriesNr).Path)
			Job.Nifti(SubjNr,SeriesNr).Path   = '';
			Job.Nifti(SubjNr,SeriesNr).Files  = {};
			Job.Nifti(SubjNr,SeriesNr).T1Text = '';
		end
	end

end

RBTM = cell(NrSubj, NrSeries);
if ~strcmp(Job.ArtDetMenu.Str{Job.ArtDetMenu.Val}, 'none')
	
	% Perform realignment and automatic artefact-detection
	NiList = makeflist(Job.Nifti);				% Update the NrSubj x NrSeries list
	WCrit  = str2num(Job.PATCHText.Str);
	if Job.ParallelBox.Val
		
		% Compute the rigid-body-transformation parameters of the DW images
		fprintf('\n==> Realigning DW images (%s)\n', datestr(now))
		[TimReq MemReq] = maxreq(NDir * 120, 0.5 * 1024^3);
		fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
		RBTM = qsubcellfun(@dd_basicproc_realignpar, stackjob(Job,'Nifti','SubjSeries'), ...
							repmat(num2cell(1:NrSubj)',1,NrSeries), repmat(num2cell(1:NrSeries),NrSubj,1), ...
							repmat(LogNames',1,NrSeries), repmat({Constr},NrSubj,NrSeries), repmat({Tuning},NrSubj,NrSeries), ...
							'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);
		[NiList{cellfun(@isempty,RBTM)}] = deal({});
		
		% Start the automatic artefact detection
		fprintf('\n==> Detecting artefacts in DW images (%s)\n', datestr(now))
		[TimReq MemReq] = maxreq(NDir*60, NDir*120*1024^2);
		fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
		ArtDFiles = qsubcellfun(@dd_patch, NiList, cell(NrSubj,NrSeries), repmat({WCrit},NrSubj,NrSeries), ...
								num2cell(true(NrSubj,NrSeries)), RBTM, repmat(LogNames',1,NrSeries), ...
								repmat(num2cell(1:NrSeries),NrSubj,1), 'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);
		
	else
		
		fprintf('\n==> Performing realignment and artefact detection (%s)\n', datestr(now))
		myset(HG, 'Visible', 'On')
		myset(HI, 'Visible', 'On')
		mywaitbar(0, HW, ['Automatic artefact detection (' num2str(NrSubj*NrSeries) ' datasets)...'])
		for SubjNr = 1:NrSubj
			for SeriesNr = 1:NrSeries
				% Compute the rigid-body-transformation parameters of the DW images
				RBTM{SubjNr,SeriesNr} = dd_basicproc_realignpar(Job, SubjNr, SeriesNr, LogNames{SubjNr}, Constr, Tuning, HWS);

				% Start the automatic artefact detection
				ArtDFiles{SubjNr,SeriesNr} = dd_patch(NiList{SubjNr,SeriesNr}, [], WCrit, true, RBTM{NrSubj,NrSeries}, LogNames{SubjNr}, SeriesNr);

				mywaitbar(((SubjNr-1)*NrSeries+SeriesNr)/numel(Job.Nifti), HW)
			end
		end
		
	end
	
	% Update the Job info
	for n = 1:NrSubj*NrSeries
		Job.Nifti(n).Files = makefcell(ArtDFiles{n});
	end

    % Now let the user interact (let's do this on his own workstation)
	if strcmp(Job.ArtDetMenu.Str{Job.ArtDetMenu.Val}, 'interactive')
		myset(HG, 'Visible', 'On')
		myset(HI, 'Visible', 'On')
		lastwarn('')				% Avoid recasting old "the job results are not yet available" warnings
        mywaitbar(0, HW, ['Interactive artefact detection (' num2str(numel(Job.Nifti)) ' subjects)...'])
		for SubjNr = 1:NrSubj
			for SeriesNr = 1:NrSeries
				
				if isempty(Job.Nifti(SubjNr,SeriesNr).Files)
					fprintf('Failed: %s\n', Job.Nifti(SubjNr,SeriesNr).Path)
					continue
				end
				
				% Start the interactive artefact detection
				fprintf('\n-> %s\n', Job.Nifti(SubjNr,SeriesNr).Path)
				delete([Job.Nifti(SubjNr,SeriesNr).Path filesep 'artc_*.*'])	% Delete results from automatic detection
				load([Job.Nifti(SubjNr,SeriesNr).Path filesep 'PATCH.mat'])
				ArtDFiles = dd_patch(PATCH, LogNames{SubjNr}, SeriesNr);
				Job.Nifti(SubjNr,SeriesNr).Files = makefcell(ArtDFiles);
				mywaitbar(((SubjNr-1)*NrSeries+SeriesNr)/numel(Job.Nifti), HW)
			end
		end
		clearvars PATCH
	end
	
else
	
	RBTM = cell(NrSubj,NrSeries);
	
end


%% --- Processing: Stage 3 (rNifti Job) ---

if Job.ParallelBox.Val
	
	myset(HI, 'Visible', 'Off')
	myset(HG, 'Visible', 'Off')

	% Realign the DW images. NB: Do not resubmit unless it takes really long (to avoid two jobs messing with the same data)
	fprintf('\n==> Realigning DW + b0 images (%s)\n', datestr(now))
	[TimReq MemReq] = maxreq(NDir*120 * NrSeries, 1024^3);
	fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, 3 * MemReq/1024^2)
	[RBTM APar D2TM] = qsubcellfun(@dd_basicproc_realign, stackjob(Job,'Nifti','Subj'), ...
									num2cell(1:NrSubj), repack(RBTM), LogNames, ...
									repmat({Constr},1,NrSubj), repmat({Tuning},1,NrSubj), ...
									'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);
	[Job.Nifti(cellfun(@isempty, RBTM),:).Files] = deal({});

	% Reslice/unwarp the DW images
	fprintf('\n==> Unwarping / resampling images (%s)\n', datestr(now))
	[TimReq MemReq] = maxreq(NDir*120 * NrSeries + 2*60^2, 1.5*1024^3);
	fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, 3 * MemReq/1024^2)
	[TgtImgs WImgs] = qsubcellfun(@dd_basicproc_warp, stackjob(Job,'Nifti','Subj'), ...
								  num2cell(1:NrSubj), LogNames, APar, D2TM, ...
								  'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);
	
	% Store the resliced/unwarped output images in the job-file (advanced use, could be handy for later use)
	for SubjNr = 1:NrSubj
		try
		if any(cellfun(@isempty,TgtImgs(SubjNr)))
			[Job.Nifti(SubjNr,:).Files] = deal({});
		end
		for SeriesNr = 1:NrSeries
			Job.Output(SubjNr,SeriesNr).TgtImgs = TgtImgs{SubjNr}{SeriesNr};
			Job.Output(SubjNr,SeriesNr).WImgs   = WImgs{SubjNr}{SeriesNr};
		end
		end
	end
	
	% Export the data for use in Camino or FSL
	fprintf('\n==> Exporting imaging data (%s)\n', datestr(now))
	[TimReq MemReq] = maxreq(NDir*10, NDir*10*1024^2);
	fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
	OK = qsubcellfun(@dd_basicproc_exportdwi, stackjob(Job,'Nifti','SubjSeries'), ...
					  repmat(num2cell(1:NrSubj)',1,NrSeries), repmat(num2cell(1:NrSeries),NrSubj,1), ...
					  'TimReq',TimReq, 'MemReq',MemReq, 'Stack',10, 'options','-l gres=bandwidth:999', 'StopOnError',false);
	
	% Estimate the diffusion tensor and its invariants
	fprintf('\n==> Estimating diffusions tensors and (normalized) invariants (%s)\n', datestr(now))
	[TimReq MemReq] = maxreq(NDir*60, NDir*20*1024^2);
	fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
	TDImgs = qsubcellfun(@dd_basicproc_estdtensor, stackjob(Job,'Nifti','SubjSeries'), ...
						 repmat(num2cell(1:NrSubj)',1,NrSeries), repmat(num2cell(1:NrSeries),NrSubj,1), ...
						 repmat(LogNames',1,NrSeries), 'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);

	% Hyper-align & normalize the tensor invariants
	fprintf('\n==> Normalizing diffusions tensor invariants (%s)\n', datestr(now))
	[TimReq MemReq] = maxreq(NDir*60 * NrSeries, NDir*20*1024^2);
	fprintf('Estimated time/job required:\t%.1f minutes\nEstimated memory/job required:\t%.1f MB\n', TimReq/60, MemReq/1024^2)
	TDImgs = qsubcellfun(@dd_basicproc_normtdimgs, stackjob(Job,'Nifti','Subj'), ...
						 num2cell(1:NrSubj), LogNames, repack(TDImgs), ...
						 'TimReq',TimReq, 'MemReq',MemReq, 'StopOnError',false);
	TDImgs = cat(1,TDImgs{:});		% Unpack the rowcell-in-cell. NB: TDImgs is an output argument
	
else
	
	myset(HI, 'Visible', 'On')
	myset(HG, 'Visible', 'On')
	
	mywaitbar(0, HW, ['Realignment, unwarping and tensor estimation (' num2str(numel(Job.Nifti)) ' subjects) ...'])
	TDImgs = cell(NrSubj, NrSeries);
	for SubjNr = 1:NrSubj
		
		% Realign the DW images
		[RBTM_ APar D2TM] = dd_basicproc_realign(Job, SubjNr, RBTM(SubjNr,:), LogNames{SubjNr}, Constr, Tuning, HWS);
		
		% Reslice/unwarp the DW images
		[TgtImgs WImgs] = dd_basicproc_warp(Job, SubjNr, LogNames{SubjNr}, APar, D2TM);

		for SeriesNr = 1:NrSeries
			
			% Store the resliced/unwarped output images in the job-file (advanced use, could be handy for later use)
			Job.Output(SubjNr,SeriesNr).TgtImgs = TgtImgs{SeriesNr};
			Job.Output(SubjNr,SeriesNr).WImgs   = WImgs{SeriesNr};
		
			% Export the data for use in Camino or FSL
			OK{SubjNr,SeriesNr} = dd_basicproc_exportdwi(Job, SubjNr, SeriesNr);
			
			% Estimate the diffusion tensor and its (normalized) invariants
			TDImgs{SubjNr,SeriesNr} = dd_basicproc_estdtensor(Job, SubjNr, SeriesNr, LogNames{SubjNr});
			
			mywaitbar(((SubjNr-1)*NrSeries+SeriesNr)/numel(Job.Nifti), HW)

		end
		
		% Hyper-align & normalize the tensor invariants. NB: TDImgs is an output argument
		TDImgs(SubjNr,:) = dd_basicproc_normtdimgs(Job, SubjNr, LogNames{SubjNr}, TDImgs(SubjNr,:));
		
	end
end
for n = 1:NrSubj*NrSeries
	if isempty(OK{n}) || ~OK{n}
		fprintf('Failed: %s\n', Job.Nifti(n).Path)
		Job.Nifti(n).Path	= '';
		Job.Nifti(n).Files	= {};
		Job.Nifti(n).T1Text = '';
	end
end

% Clean-up the intermediate data files
Job = dd_basicproc_cleanup(Job);
myclose(HI, HG, HW, HWS)


%% =-=-=-=-=-=-=-=-=-=- Auxillary Functions -=-=-=-=-=-=-=-=-=-=-=

function dd_basicproc_bvalvec2mat(Job, n)

for m = 1:numel(Job.Nifti(n).Files)
	DWIFiles{m} = fullfile(Job.Nifti(n).Path, Job.Nifti(n).Files{m});
end
try
	dti_get_dtidata(DWIFiles{1});
catch
	warning('No standard DWI information found in %s. Will try to use bval and bvec files', DWIFiles{1})
	if spm_select('List', Job.Nifti(n).Path, 'bvec\.txt|bvec$|bvecs$')
		dd_bvalvec2mat(DWIFiles);
	else
		error('No DWI information found in %s', Job.Nifti(n).Path)
	end
end


function Subj = repack(SubjSeries)

for SubjNr = 1:size(SubjSeries,1)
	SubjSeries{SubjNr,1} = SubjSeries(SubjNr,:);	% Create a cellrow-in-cell
end
Subj = SubjSeries(:,1)';


function QJob = stackjob(Job, Mode, Stacking)

% Make a stacked Job cell array (for qsubcellfun) with minimal memory footprint

switch Mode
	case 'DICOM'
		Job = rmfield(Job, 'Nifti');
	case 'Nifti'
		Job = rmfield(Job, 'DICOM');
end

% Stack the Job structure (memory footprint ~ NrSubj^2)
switch Stacking
	case 'Subj'
		QJob = repmat({Job},1,size(Job.(Mode),1));	% Stacking per Subject (for processing entire Series-sets)
	case 'SubjSeries'
		QJob = repmat({Job},size(Job.(Mode)));		% Stacking per Subject x Series
	case 'LinInd'
		QJob = repmat({Job},1,numel(Job.(Mode)));	% Stacking per Subject, Series
end

% Remove all redundant non-central/off-diagonal data
for n = 1:numel(QJob)
	for m = 1:numel(QJob)
		if n~=m
			switch Stacking
				case {'SubjSeries','LinInd'}		% Remove all non-central cells
					QJob{n}.(Mode)(m).Path   = '';
					QJob{n}.(Mode)(m).Files  = {};
					QJob{n}.(Mode)(m).T1Text = '';
					if iscell(Job.BETOpts.Str)		% Use individual BET settings
						QJob{n}.BETOpts.Str{m} = {};
					end
					if isfield(QJob{n}, 'Output')
						QJob{n}.Output(m).TgtImgs = {};
						QJob{n}.Output(m).WImgs	  = {};
					end
				case  'Subj'						% Remove all non-central rows
					[QJob{n}.(Mode)(m,:).Path]   = deal('');
					[QJob{n}.(Mode)(m,:).Files]  = deal({});
					[QJob{n}.(Mode)(m,:).T1Text] = deal('');
					if iscell(Job.BETOpts.Str)		% Use individual BET settings
						[QJob{n}.BETOpts.Str{m,:}] = deal({});
					end
					if isfield(QJob{n}, 'Output')
						[QJob{n}.Output(m,:).TgtImgs] = deal({});
						[QJob{n}.Output(m,:).WImgs]	  = deal({});
					end
			end
		end
	end
end


function mywaitbar(varargin)

% Robust against closed figures
if ishandle(varargin{2})
    waitbar(varargin{:})
end


function myset(varargin)

% Robust against closed figures
if ishandle(varargin{1})
	set(varargin{:})
end


function OK = tryset(handles, Job, Item)

% Robust against loading incompatible fields from old Job-files
try
	if strfind(Item, 'Menu')
		NewVal = find(strcmpi(get(handles.(Item),'String'), Job.(Item).Str(Job.(Item).Val)));
		if isempty(NewVal), error, end
		set(handles.(Item), 'Value',NewVal)
	elseif strfind(Item, 'Box')
		set(handles.(Item), 'Value',Job.(Item).Val)
	elseif ~isempty(strfind(Item,'Text'))
		set(handles.(Item), 'String',Job.(Item).Str)
	else
		error('TrySet:Item', 'Unknown Item: %s', Item)
	end
catch Exception
	if strcmp(Exception.identifier, 'TrySet:Item')
		rethrow(Exception)
	else
		warning(['Could not set job-item: ' Item])
	end
end

if nargout
	if exist('Exception','var')
		OK = false;
	else
		OK = true;
	end
end


function myclose(varargin)

% Robust against closed figures
for n = 1:numel(varargin)
	if ishandle(varargin{n})
		close(varargin{n})
	end
end


function update_data_panel(handles, SubjNr, SeriesNr)

if nargin<3
    SeriesNr = get(handles.MeasurementMenu, 'Value');
end
if nargin<2 || isempty(SubjNr)
    SubjNr = get(handles.SubjMenu, 'Value');
end

% Get datatype
DataTypes = get(handles.DataMenu, 'String');
Selected  = get(handles.DataMenu, 'Value');

% Update data set
DirText  = [];
FileText = [];
T1Text   = '';
[NrSubj NrSeries] = size(handles.(DataTypes{Selected}));
if NrSubj
	if isinf(SubjNr)		% Take the last subject
		SubjNr = NrSubj;
	end
	if isinf(SeriesNr)		% Take the first series
		SeriesNr = 1;
	end
	DirText  = handles.(DataTypes{Selected})(SubjNr, SeriesNr).Path;
	FileText = handles.(DataTypes{Selected})(SubjNr, SeriesNr).Files;
	T1Text   = handles.(DataTypes{Selected})(SubjNr, SeriesNr).T1Text;
	if ~isempty(T1Text) && size(T1Text,1)>1
		T1Text = [deblank(T1Text(1,:)) ',*'];
	end
end
set(handles.DirText,  'String', DirText)
set(handles.FileText, 'String', FileText, 'Value', 1)
set(handles.T1Text,   'String', T1Text)

% Update subject number
Subj = '-'; Subj(1) = [];		% size 1-by-0
for n = 1:NrSubj
    Subj{n} = num2str(n);
end
if isinf(SubjNr)
    SubjNr = 1;
end
set(handles.SubjMenu, 'String', Subj, 'Value', SubjNr)

% Update measurement number
for n = 1:max(1,NrSeries)
    Series{n} = num2str(n);
end
if isinf(SeriesNr)
	SeriesNr = 1;
end
set(handles.MeasurementMenu, 'String', Series, 'Value', SeriesNr)


function update_longitudinal_panel(handles)

if numel(get(handles.MeasurementMenu, 'String'))>1
	if get(handles.T1Box,'Value')
		set(handles.T1RealignBox, 'Visible', 'on')
		if get(handles.PEUnwarpBox,'Value')
			set(handles.SinglePEUnwarpBox, 'Visible', 'on')
		else
			set(handles.SinglePEUnwarpBox, 'Visible', 'off')
		end
	else
		set(handles.T1RealignBox, 'Visible', 'off')
		set(handles.SinglePEUnwarpBox, 'Visible', 'off')
	end
	set(handles.HyperalignMenu, 'Visible', 'on')
	set(handles.HyperalignText, 'Visible', 'on')
else
	set(handles.T1RealignBox, 'Visible', 'off')
	set(handles.SinglePEUnwarpBox, 'Visible', 'off')
	set(handles.HyperalignMenu, 'Visible', 'off')
	set(handles.HyperalignText, 'Visible', 'off')
end


function update_preprocessing_panel(handles)

SeriesNr = get(handles.MeasurementMenu, 'Value');
SubjNr	 = get(handles.SubjMenu, 'Value');
if ~isempty(handles.BETOpts_.Str)
	set(handles.BETOpts, 'String', handles.BETOpts_.Str{SubjNr,SeriesNr})
	set(handles.BETMenu, 'Value',  handles.BETMenu_.Val(SubjNr,SeriesNr))
end

if get(handles.T1Box,'Value')
	set(handles.T1Text, 'Visible', 'on')
	set(handles.PEUnwarpBox, 'Visible', 'on')
	if get(handles.PEUnwarpBox, 'Value')
		set(handles.PEUnwarpText, 'Visible', 'on')
	else
		set(handles.PEUnwarpText, 'Visible', 'off')
	end
else
	set(handles.T1Text, 'Visible', 'off')
	set(handles.PEUnwarpBox, 'Visible', 'off')
	set(handles.PEUnwarpText, 'Visible', 'off')
end

Method = cellstr(get(handles.DenoisingMenu,'String'));
if strcmp(Method{get(handles.DenoisingMenu,'Value')}, 'none')
	set(handles.RicianBox, 'Visible', 'off')
else
	set(handles.RicianBox, 'Visible', 'on')
end

% If 'none' is chosen, check if PATCH estimation is used
Options = cellstr(get(handles.ArtDetMenu, 'String'));		% Options as cell array
if strcmp(Options{get(handles.ArtDetMenu, 'Value')}, 'none')
    set(handles.PATCHText, 'Visible', 'off')
	EstMethod = cellstr(get(handles.EstMenu,'String'));		% EstMenu contents as cell array
	if strcmp(EstMethod{get(handles.EstMenu,'Value')}, 'PATCH')
		warndlg('Estimation method set to ''GLM''')
		set(handles.EstMenu, 'Value', find(strcmp('GLM', get(handles.EstMenu, 'String'))))	% Set method to GLM
	end
else
    set(handles.PATCHText, 'Visible', 'on')
end

if get(handles.PEUnwarpBox,'Value') && get(handles.T1Box,'Value')
	set(handles.PEUnwarpText, 'Visible', 'on')
else
	set(handles.PEUnwarpText, 'Visible', 'off')
end


function update_tensor_panel(handles)

Options = cellstr(get(handles.EstMenu, 'String'));				% Options as cell array
if strcmp(Options{get(handles.EstMenu, 'Value')}, 'PATCH')
	ArtDetMethod = cellstr(get(handles.ArtDetMenu,'String'));	% ArtDetMenu contents as cell array
	if strcmp(ArtDetMethod{get(handles.ArtDetMenu,'Value')}, 'none')
		warndlg('Artefact detection method set to ''automatic''')
		set(handles.ArtDetMenu, 'Value', find(strcmp('automatic', get(handles.ArtDetMenu, 'String'))))
	end
	set(handles.PATCHText, 'Visible', 'on')
end
if strcmp(Options{get(handles.EstMenu, 'Value')}, 'none')
	set(handles.MDBox, 'Visible', 'off')
	set(handles.NormDBox, 'Visible', 'off')
	set(handles.NormABox, 'Visible', 'off')
	set(handles.FABox, 'Visible', 'off')
	set(handles.RABox, 'Visible', 'off')
	set(handles.ModeBox, 'Visible', 'off')
	set(handles.EigBox, 'Visible', 'off')
	set(handles.RDBox, 'Visible', 'off')
	set(handles.DottedLine, 'Visible', 'off')
	set(handles.NormBox, 'Visible', 'off')
else
	set(handles.MDBox, 'Visible', 'on')
	set(handles.NormDBox, 'Visible', 'on')
	set(handles.NormABox, 'Visible', 'on')
	set(handles.FABox, 'Visible', 'on')
	set(handles.RABox, 'Visible', 'on')
	set(handles.ModeBox, 'Visible', 'on')
	set(handles.EigBox, 'Visible', 'on')
	set(handles.RDBox, 'Visible', 'on')
	set(handles.DottedLine, 'Visible', 'on')
	set(handles.NormBox, 'Visible', 'on')
end


function [handles SubjNr] = appendjob(handles, JobFile)

% Get existing data type and size
DataTypes		  = get(handles.DataMenu, 'String');
Selected		  = get(handles.DataMenu, 'Value');
[Subjects Series] = size(handles.(DataTypes{Selected}));

% Load the new data and check if there is data and if the nr of measurements is the same
load(JobFile)
[NewSubjects NewSeries] = size(Job.(DataTypes{Selected}));
if NewSubjects==0
	warndlg(['No ' DataTypes{Selected} ' files found in Job-file'], ...
			'Append data sets in Job')
	return
end
if (NewSeries~=Series) && Subjects>0
	warndlg(['Number of measurements / time points in Job-file (' num2str(NewSeries) ...
			 ') does not match the existing number of measurements / time points (' num2str(Series) ')'], ...
			'Append data sets in Job')
	SubjNr = 1;
	return
end

% Append the new DICOM, Nifti and BET data
if isempty(handles.DICOM)
	handles.DICOM = Job.DICOM;
elseif size(Job.DICOM,1) == NewSubjects		% If the Job contains data add it
	handles.DICOM(Subjects + (1:NewSubjects), 1:NewSeries) = Job.DICOM;
elseif isstruct(handles.DICOM)				% Ooops, I forgot what this is for?
	handles.DICOM(Subjects + (1:NewSubjects), 1:NewSeries) = struct('Path',[], 'Files',[], 'T1Text',[]);
end
if isempty(handles.Nifti)
	handles.Nifti = Job.Nifti;
elseif size(Job.Nifti,1) == NewSubjects		% If the Job contains data add it
	handles.Nifti(Subjects + (1:NewSubjects), 1:NewSeries) = Job.Nifti;
elseif isstruct(handles.Nifti)				% Ooops, I forgot what this is for?
	handles.Nifti(Subjects + (1:NewSubjects), 1:NewSeries) = struct('Path',[], 'Files',[], 'T1Text',[]);
end
handles.BETMenu_.Val(Subjects + (1:NewSubjects), 1:NewSeries) = Job.BETMenu.Val;
handles.BETOpts_.Str(Subjects + (1:NewSubjects), 1:NewSeries) = cellstr(Job.BETOpts.Str);
set(handles.BETBox,'Value', 1)				% This may not always be necessary
SubjNr = Subjects + NewSubjects;


function FNames = makefcell(FList)
%
% Input:  filelist (SPM-style)
% Output: filecell (Job-style)

if isempty(FList)
	FNames = '';
	return
end

Path0 = fileparts(FList(1,:));
for n = 1:size(FList,1)
    [Path FName Ext] = fileparts(FList(n,:));
    FNames{n}		 = [FName strtok(Ext,',')];		% Throw away frame-info
    if ~strcmp(Path, Path0)
        error('Files from 1 series must be in the same directory')
    end
end


function CellList = makeflist(DataStruct)
%
% Input:  Job.DataStruct (.Path & .Files)
% Output: cell-in-cellarray of size Subjects*Series (for passing data to dd_patch)

CellList = cell(size(DataStruct));
for n = 1:numel(DataStruct)
	FList = cell(1,numel(DataStruct(n).Files));
	for m = 1:numel(DataStruct(n).Files)
		FList{m} = fullfile(DataStruct(n).Path, DataStruct(n).Files{m});
	end
	CellList{n} = {char(FList)};
end


function PO = prepend(PI, pre)
%
% Input: filelist (SPM-style) or cellarray

for n = 1:size(PI,1)
    [pth,nm,xt] = fileparts(char(PI(n,:)));
    PO{n}       = fullfile(pth, [pre nm xt]);
end
if ischar(PI)
	PO = char(PO);
end


function OK = check_t1(hObject, handles)

DataTypes = get(handles.DataMenu, 'String');
Selected  = get(handles.DataMenu, 'Value');

% Get the files from the user
switch DataTypes{Selected}
    case 'DICOM'
        Data = handles.DICOM(:);
    case 'Nifti'
        Data = handles.Nifti(:);
end
if numel(Data)==0 || (~isempty(Data(1).T1Text) && all(~strcmp('', {Data.T1Text}))) % OLD: isempty(strmatch('', {Data.T1Text}, 'exact')))
    OK = true;
    return    
end

switch questdlg(['Some data sets have unspecified T1-reference images. Do you want to ' ...
                'continue or do you want to specify them now?'], 'Question', 'Specify', 'Continue', 'Cancel', 'Specify');
    
    case 'Specify'

        % Check the data sets for empty T1Text-fields
        for n = 1:numel(Data)
            if isempty(Data(n).T1Text)
                switch DataTypes{Selected}
                    case 'DICOM'
                        [T1Text Sts] = spm_select(Inf, 'any', ...
                            ['Select T1 reference scan (DICOM or Nifti) for data-set ' num2str(n)], {''}, Data(n).Path,'\.nii$|\.img$|(?i)\.ima$|(?i)\.dcm$');
                        handles.DICOM(n).T1Text = T1Text;
                    case 'Nifti'
                        [T1Text Sts] = spm_select(1, 'image', ['Select T1 reference scan for data-set ' num2str(n)], {''}, Data(n).Path);
                        handles.Nifti(n).T1Text = strtok(T1Text, ',');  % Throw away frame-info;
                end
                if Sts==false, break, end
            end
        end
        OK = false;  % false = return to the gui
        
        % Update data set
		guidata(hObject, handles)
		update_data_panel(handles)

    case 'Continue'

        OK = true;
        
    otherwise
        
        OK = [];

end

function [TimReq MemReq] = maxreq(TimReq, MemReq)

% TODO: retrieve values from qstat -q
MaxTime = 48*60*60;
MaxMem  = 256*1024^3;

if TimReq > MaxTime
	warning('Required time %i is larger than the queue maximum %i', TimReq, MaxTime)
	TimReq = MaxTime;
end
if MemReq > MaxMem
	warning('Required memory %i is larger than the queue maximum %i', MemReq, MaxMem)
	MemReq = MaxMem;
end
