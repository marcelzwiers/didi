function varargout = dd_patch(varargin)

% DD_PATCH automatically detects excessive inplane residual errors (e.g. due to
% cardiac or subject motion) in DW images (b > 50 s/mm^2). Artefact corrected
% files ('artc_*') and PATCH info will be saved to disk ('PATCH.*').
%
% A full description of the method can be found in:
% Zwiers, M.P. (2010). Patching cardiac and head motion artefacts in
% diffusion-weighted images. Neuroimage doi:10.1016/j.neuroimage.2010.06.014
%
%
% Requirements:
% ~~~~~~~~~~~~~
% - SPM8/12 and FSL (optional, see commandline usage of Mask)
%
%
% GUI Usage:
% ~~~~~~~~~~
% Running DD_PATCH by itself (no input arguments) creates a new gui or
% raises the existing one. The upper panels show the minimal voxel
% (left) and sum (right) of the in plane residual error according to the
% diffusion tensor model.
%   Slices for which the residual errors exceed the threshold (dashed line;
% can be set with [Threshold]) are considered outliers and are shown in
% the lower panels. You can navigate through these outliers by using the
% slider. The goal is to decide which image is best, the uncorrected image
% (left), the partially corrected image (middle) or the slice corrected
% image (right). In addition to browsing just through the outliers, you can
% also browse through the slices ([z+]/[z-]), through the different
% diffusion directions ([b+]/[b-]), or directly pick a particular slice
% ([set]). Note that navigation functionality is also available by using
% your keyboard: <Left Arrow>/<Right Arrow> for browsing through the
% detected outliers, <Arrow Up>/<Arrow Down> for browsing in the
% z-direction, and the <Ctrl>+<Left/Right Arrow> for browsing through
% neighbouring DW acquisitions (diffusion directions).
%   When you are satisfied with all of the image selections you can press
% [Accept]. If you find the algorithm not sensitive enough or too
% sensitive, or you feel that the data is corrupted by vibrations from the
% x-gradient of the scanner, you can restart the outlier detection process
% ([Recalculate]) using different weighting criteria (see also the 'WCrit'
% description below).
%
%
% Commandline Usage:
% ~~~~~~~~~~~~~~~~~~
%
% FUNCTION [DTIFiles PATCH] = DD_PATCH(DTIFiles, Mask, WCrit, Auto, AffM)
%
% INPUT
%   DTIFiles - Cell containg char array of DTI filenames (optional)
%              NB: Matlab crashes if it's a char array
%   Mask     - Mask or name of Mask-file of voxels to include. If left
%              empty an automatically generated mask is used (requires FSL bet)
%   WCrit    - Scaling vector with 5 entries:
%              WCrit(1) = Voxel (cardiac) weighting. This value ranges
%                         from zero to infinity: zero gives equal weighting
%                         and larger values give less weighting to larger
%                         residuals (e.g. 0.3 supresses residuals larger
%                         than 3 times the standard deviation). Negative
%                         values evoke the use of image closure of the
%                         weights instead of the normal median filtering of
%                         the residuals (see paper for details).
%              WCrit(2) = Slicewise (head motion) weighting. This value
%                         ranges in the same way as WCrit(1)
%              WCrit(3) = Logarithmic intensity weighting (Salvador et al., 
%                         2005). The weight ranges from 0 (no weighting) to
%                         1 (full weighting)
%              WCrit(4) = Weigting value for x-vibration affected volumes.
%                         The weight ranges from 0 to 1
%              WCrit(5) = Threshold for applying x-vibration weighting. The
%                         value represents the x-component of the diffusion
%                         gradient and ranges from 0 to 1. A value of zeros
%                         includes all directions and a value of one
%                         includes only pure x-directions.
%              When WCrit = [] the default is used: [0.15 0.15 0 1 0.7]
%              When WCrit = NaN the GMM method is used (Basser et al., 2006)
%   Auto     - Skips user interaction if true (default false)
%   AffM     - Rotation matrices used for realignment
%
% OUTPUT
%   DTIFiles - char array of artefact-corrected DTI filenames
%   PATCH    - Struct with fields containing information of IRLS regression
%              with fields: 'input' (DTIFiles), 'MaskFile' (Mask), 'b0File',
%			   'WCrit', 'WVox' (voxel weightings), 'WSlc' (slice-wise weightings),
%              'Thresh' (handles.Thresh), 'Outliers' (user selections)
%
% 
% FUNCTION [..] = DD_PATCH(PATCH)
%
% Alternative use of DD_PATCH in which the calculation restarts from the
% endpoint of a previous calculation
%
% INPUT
%   PATCH    - IRLS struct from previous DD_PATCH run (see also above)
%
% _________________________________________________________________________
% DD_PATCH was developed by Marcel Zwiers. Please report all bugs to:
% m.zwiers@donders.ru.nl
%
% Ver 4.0, 23/04/2010. For GPLv3 license see: www.ru.nl/neuroimaging/diffusion
%
% See also: DD_BASICPROC.

% DD_PATCH M-file for dd_patch.fig
%      DD_PATCH, by itself, creates a new DD_PATCH or raises the existing
%      singleton*.
%
%      H = DD_PATCH returns the handle to a new DD_PATCH or the handle to
%      the existing singleton*.
%
%      DD_PATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DD_PATCH.M with the given input arguments.
%
%      DD_PATCH('Property','Value',...) creates a new DD_PATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dd_patch_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dd_patch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dd_patch

% Last Modified by GUIDE v2.5 05-Aug-2009 15:32:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dd_patch_OpeningFcn, ...
                   'gui_OutputFcn',  @dd_patch_OutputFcn, ...
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


%% --- Executes just before dd_patch is made visible.
function dd_patch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dd_patch (see VARARGIN)
% 
% USAGE I:
% dd_patch_OpeningFcn(hObject, eventdata, handles)
%	-> all input arguments are acquired interactively
%
% USAGE II:
% dd_patch_OpeningFcn([..], PATCH, [LogName], [SeriesNr])
%	-> Used for pre-computed interactive sessions
%	-> Info is copied from previously saved PATCH file
%	-> LogName is used in spm_print (handy for distributed computing)
%	-> SeriesNr is used in spm_print (handy for distributed computing)
%
% USAGE III:
% dd_patch_OpeningFcn([..], DTIFiles, Mask, WCrit, Auto, AffM, [LogName], [SeriesNr])
%	-> Default commandline use
%   -> Mask is a string (=> dd_basicproc) or matrix.
%	-> LogName is used in spm_print (handy for distributed computing)
%	-> SeriesNr is used in spm_print (handy for distributed computing)
%
% USAGE IV:
% dd_patch_OpeningFcn([..], DTIFiles, Mask, WCrit, Auto, AffM, Outliers)
%	-> Used by GUI => Recalc
%	-> Info from previous call is passed via handles.[*] (no need to load
%	   or compute them again = faster).
%   -> Mask is a matrix

% Add the path to the Diffusion toolbox
AddPath = false;
WhatdTB = what('Diffusion');     % Check if the Diffusion-toolbox on the Matlab-path
if numel(WhatdTB)==0
    xTB = spm('TBs');            % Check if the Diffusion-toolbox is installed in SPM
    dTB = find(strcmp('Diffusion', {xTB.name}));
    addpath(genpath(xTB(dTB).dir))
    AddPath = true;
end

% -- Initialize some variables
WVox	= [];
WSlc	= [];
B0		= 50;											% Threshold for selecting b0 images
Thresh	= [0.5 0.5];									% Default values
WCrit	= [0.15 0.15 0 1 0.7];							% Default values
tzyx	= [4 3 2 1];									% This dim-order is convenient for later calculations
Keep	= false;
R		= 2;											% Radius of the Mask dilation kernel
Kernel  = ones(2*R+1, 2*R+1, 2*R+1);
[I J K] = ind2sub(size(Kernel), find(Kernel));
Kernel(sqrt((I-R-1).^2 + (J-R-1).^2 + (K-R-1).^2) > R) = 0;

% -- Parse the input variables
IsRecalc  = false;
[S LWarn] = mywarning('Off', 'SPM:noDisplay');			% Else spm_figure gives warning when using qsubcellfun
if (numel(varargin)>5) && isnumeric(varargin{6})		% Usage IV (= Recalculation from the GUI)
	IsRecalc = true;
elseif (numel(varargin)>5) && ischar(varargin{6}) && ~isempty(varargin{6})	% Usage III + distributed computing
	spm_jobman('initcfg');
	HG = spm_figure('GetWin', 'Graphics');
	setappdata(HG, 'LogName', varargin{6});
	if numel(varargin)>6
		setappdata(HG, 'SeriesNr', varargin{7});
	end
end
if numel(varargin)>=1 && numel(varargin)<=3 && isfield(varargin{1}, 'WVox')	% Usage II
    PATCH = varargin{1};
	HG	  = spm_figure('FindWin', 'Graphics');
	if numel(varargin)>1 && ishandle(HG)
		setappdata(HG, 'LogName', varargin{2});
		if numel(varargin)>2
			setappdata(HG, 'SeriesNr', varargin{3});
		end
	end
elseif numel(varargin)==0								% Usage I => get PATCH or DTIFiles
	set(spm_figure('GetWin','Graphics'), 'visible','off')
	set(spm_figure('GetWin','Interactive'), 'visible','off')
	ValidInput = false;
	while ~ValidInput
		SelFiles = spm_select(Inf, 'any', 'Select the DWI files or PATCH.mat from a previous session', ...
							  {''}, pwd, 'PATCH.mat|\.img$|\.nii$');
		if size(SelFiles,1)==0
			return
		elseif size(SelFiles,1)==1
			if regexp(SelFiles, 'PATCH.mat\>')			% An PATCH-file was selected
				load(SelFiles)							% Load the PATCH-file
				ValidInput = true;
				Keep	   = true;						% Keep selected outliers (=> UsrOutliers)
			else
				uiwait(warndlg('More than 6 DWI volumes required'))
			end
		elseif size(SelFiles,1)>6
			varargin{1} = SelFiles;						% DWI-files were selected
			ValidInput  = true;
		else
			uiwait(warndlg('Select more than 6 DWI files or 1 PATCH.mat file'))
		end
	end
end
mywarning(S, LWarn)
if exist('PATCH', 'var')								% Copy info from the PATCH variable (Usage I or II)
    DTIFiles 	= PATCH.input;							% Load the DTI volumes from disk
    Mask		= PATCH.MaskFile;						% Load the Mask volume from disk
	b0File		= PATCH.b0File;
    WVox		= PATCH.WVox;
    WSlc		= PATCH.WSlc;
	WLog		= PATCH.WLog;
    Thresh		= PATCH.Thresh;
    UsrOutliers	= PATCH.Outliers;
	for n = 1:numel(PATCH.WCrit)
		WCrit(n) = PATCH.WCrit(n);
	end
	if isfield(PATCH, 'AffM')
		AffM = PATCH.AffM;
	else
		AffM = [];
	end
else													% Get the info from the input
	if isempty(varargin{1})
		handles.close = true;							% Used in dd_patch_outputfcn
		guidata(hObject, handles);
		return
	end
	DTIFiles = char(varargin{1});
	if numel(varargin)<2
		Mask = spm_select([0 1], 'image', 'Select a Mask file (leave empty for automaticly generated mask)', {''}, pwd,  'mask.*');
	else
		Mask = varargin{2};
	end
	if numel(varargin)>=3
		for n = 1:numel(varargin{3})
			WCrit(n) = varargin{3}(n);
		end
	end
    if numel(varargin)>=4 && varargin{4}				% Auto = No user-interaction
        handles.close = true;							% Used in dd_patch_outputfcn
    end
	if numel(varargin)<5 || isempty(varargin{5})
		AffM = [];
	else
		AffM = varargin{5};
	end
	WLog = [];
end

% -- Get necessary variables
fprintf('\n%s\nRunning ''PATCH''\n', repmat('-',[1 40]))
if IsRecalc												% Usage IV (Previous info passed via handles (no need to load them again=faster))
	if any(varargin{6}>=2)								% 1=DWI, 2=PCDWI, 3=CDWI
		UsrOutliers = varargin{6};
		Keep = strcmp(questdlg('Do you want to keep the current selection?', ...
							   'Recalculate residuals', 'Yes', 'No', 'No'), 'Yes');
	end
	tzyx     = handles.tzyx;
	DWISel   = handles.DWISel;
	q        = handles.q;
	Vol      = handles.Vol;
	MaskFile = handles.MaskFile;
	b0File	 = handles.b0File;
	BBox     = handles.BBox;
	WVox     = handles.WVox;
	WSlc     = handles.WSlc;
	WLog	 = handles.WLog;
	Thresh   = handles.Thresh;
	VolHdr	 = spm_vol(DTIFiles(DWISel,:));
else													% Load the info (NB: mask avoids log of zero)
    disp('Loading data...')
	for n = 1:size(DTIFiles,1)
		D(n) = dti_get_dtidata(DTIFiles(n,:));
	end
    DWISel = [D(:).b] > B0;								% We're not interested in the b0-volumes
    q      = vertcat(D(DWISel).g);						% Matrix containing DTI gradient directions
	VolHdr = spm_vol(DTIFiles(DWISel,:));
	try													% Assume all orientations are identical
		Vol = permute(spm_read_vols(VolHdr), tzyx);
	catch												% Try reading them otherwise
		disp('WARNING: Ignoring orientation mismatches and proceeding...')
		Vol = zeros([numel(VolHdr) VolHdr(1).dim(3:-1:1)]);
		for n = 1:numel(VolHdr)
			Vol(n,:,:,:) = permute(spm_read_vols(VolHdr(n)), tzyx(2:end));
		end
	end
	if isempty(Mask)
		[Mask b0File] = dd_getmeanmask(DTIFiles(~DWISel,:), true);
	end
	if ischar(Mask)
		MaskFile = Mask;
		Mask     = spm_read_vols(spm_vol(Mask));
	end
	Mask	= permute(Mask>0.5, tzyx(2:end));			% Use same dimorder as Vol
    Myi     = find(any(any(Mask,3)));
    Mxi     = find(any(any(Mask,2)));
    BBox    = [min(Myi) max(Myi) min(Mxi) max(Mxi)];	% Bounding box Mask (used for display)
    BBox    = BBox + round(0.06*max(BBox))*[-1 1 -1 1];	% Add a border
    BBox(BBox < 1) = 1;									% Clip to avoid errors
    BBox(2) = min([BBox(2) size(Mask,2)]);				% Clip to avoid errors
    BBox(4) = min([BBox(4) size(Mask,3)]);				% Clip to avoid errors
end
if ~exist('b0File','var') || isempty(b0File)			% Make sure we have a b0 image
	b0File = DTIFiles(~DWISel,:);
	b0File = b0File(1,:);
end
b0Vol = permute(spm_dilate(spm_read_vols(spm_vol(b0File)), Kernel), tzyx(2:end));

% WCrit(1) is a little bit misused here (may be improved in future version)
if isnan(WCrit(1))
    Method = 'GMM';
	Thresh = sqrt([0.1 0.1]);  % Three times the SD (NB WVox must be scaled to 1)
else
    Method = 'PATCH';
end
PathStr = fileparts(VolHdr(1).fname);
if length(PathStr)>60, PathStr = ['...' PathStr(end-57:end)]; end
set(hObject, 'Name', ['DTI artefact detection 4.0 (' Method '): ' PathStr], 'NumberTitle', 'Off')

% Be nice to our RAM
Vol  = single(Vol);
WVox = single(WVox);
WLog = single(WLog);

% -- Compute the WLS residuals
disp('Analyzing residuals...')
[handles.CVol WVox WSlc WLog handles.SclF] = ...
    runpatch(Vol, VolHdr(1), Mask, b0Vol, q, WCrit, WVox, WSlc, WLog, Method, IsRecalc, AffM, DTIFiles(DWISel,:));
handles.CVol(:,~Mask) = Vol(:,~Mask);		% Restore values outside the mask
handles.PCMask = WVox < Thresh(1);
WMin = min(WVox(:,:,:), [], 3);				% Used for detection

% -- Detect the outliers
Outliers = [WMin(:)<Thresh(1) WSlc(:)<Thresh(2)];
if Keep
    handles.Outliers = any(Outliers,2)*1;	% 1=DWI (Let the user decide to use the new outliers)
    handles.Outliers(UsrOutliers>1) = UsrOutliers(UsrOutliers>1); % Restore the user values
else
    handles.Outliers = Outliers(:,1)*2;		% 2=PCDWI
    handles.Outliers(Outliers(:,2)) = 3;	% 3=CDWI (overrule PCDWI)
end
NrOutliers = numel(find(handles.Outliers));
fprintf('Detected voxel outliers:  %g/%g\n', sum(handles.Outliers==2), numel(handles.Outliers))
fprintf('Detected slicewise outliers: %g/%g\n', sum(handles.Outliers==3), numel(handles.Outliers))

% -- Show the results in the GUI
wplot(handles.axes_MinRes, WMin, Thresh(1), 'b', ' Voxel ')
ylabel(handles.axes_MinRes, 'Minimum weight')
wplot(handles.axes_SumRes, WSlc, Thresh(2), 'r', ' {\Sigma} Slice ')
set(handles.axes_SumRes, 'YTickLabel', '')

% -- Add userdata to handles structure (NB: some have already been assigned above)
handles.Vol       = Vol;
handles.Mask      = Mask;
handles.MaskFile  = MaskFile;
handles.b0File	  = b0File;
handles.Thresh    = Thresh;
handles.WCrit     = WCrit;
handles.WVox      = WVox;
handles.WMin      = WMin;
handles.WSlc      = WSlc;
handles.WLog	  = WLog;
handles.AffM      = AffM;
handles.BBox      = BBox;
handles.q         = q;
handles.input     = DTIFiles;
handles.DWISel    = DWISel;
handles.tzyx      = tzyx;
if NrOutliers == 0
    % Make image-panels invisible
    initimpanels(handles, 'Off')
    handles.Ind = 1;    % Initialize selected image
    set(handles.text_File, 'String', 'No outliers detected...')
end
if NrOutliers <= 1  % There is outlier slider to show
    set(handles.slider_Ind, 'Value', 1, 'Min', 1-eps, 'Max', 1, 'SliderStep', [1 1])
else
    set(handles.slider_Ind, 'Value', 1, 'Min', 1, 'Max', NrOutliers, ...
        'SliderStep', [1 1]/(NrOutliers-1))
end
guidata(hObject, handles);

% Draw first slice in the GUI
if NrOutliers>=1 && ~isfield(handles, 'close')
    slider_Ind_Callback(handles.slider_Ind, [], handles)
end

% Stop user from hitting buttons
if isfield(handles, 'close')
    set(handles.pushbutton_accept, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
    set(handles.pushbutton_recalc, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
end

% Remove the added path
if AddPath
    rmpath(genpath(xTB(dTB).dir))
end

%%


% --- Executes during object creation, after setting all properties.
function dd_patch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dd_patch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Check for screensize and licence
Units = get(hObject, 'Units');			% Store the units of the GUI
set(hObject, 'Units', get(0, 'Units'))	% Set them the same as root
GUISz = get(hObject, 'OuterPosition');	% Get the size of the GUI
set(hObject, 'Units', Units)			% Restore the units of the GUI
ScrSz = get(0, 'ScreenSize');
if ScrSz(4) < GUISz(4) + 30 + 40		% 30 = menubar, 44 = KDE-taskbar
%     warndlg(sprintf(['Your screen size may be too small (resolution too low) ' ...
% 					 'to properly display the gui-window (%g x %g)'], GUISz(3:4)), ...
% 					 'Warning DD_PATCH', 'modal')
    fprintf(['\nWarning: Your screensize may be too small (resolution too low) ' ...
			 'to properly display the gui-window (%g x %g)\n'], GUISz(3), GUISz(4))
end

% Add a menubar togglebutton to the toolbar
MIcon = repmat(shiftdim([0.92 0.91 0.84],-1), [16 16 1]);
MIcon(2:4,2:15,:)		= 0.8;	% menubar background
MIcon(3,[4:8 12:15],:)	= 0;	% menubar text
MIcon([1 5],:,:)		= 0.5;	% border
MIcon(:,[1 16],:)		= 0.5;
MIcon(6:14,[4 11],:)	= 0.5;
MIcon(15,4:11,:)		= 0.5;
MIcon(6:14,5:10,:)		= 1;	% menu
MIcon(7,5:8,:)			= 0;	% menu text
MIcon(9,5:9,:)			= 0;
MIcon(11,5:6,:)			= 0;
MIcon(13,5:7,:)			= 0;
Htb = findall(hObject,'Type','uitoolbar');
uitoggletool(Htb, 'CData', MIcon, 'TooltipString','Show menubar',...
             'HandleVisibility','callback', 'Separator','on', ...
             'OffCallBack', 'set(gcbf,''MenuBar'',''none'')', ...
             'OnCallBack', 'set(gcbf,''MenuBar'',''figure'')');

% % Change the 'Save' button to 'Save As (jpeg)' to avoid overriding dd_patch.fig
Hsv = findall(hObject, 'Tag', 'Standard.SaveFigure');
set(Hsv, 'ClickedCallback', 'filemenufcn(gcbf,''FileSaveAs'')', ...
         'TooltipString', 'Save Figure As (NB saving as dd_patch.fig is not recommended')
Types  = localexporttypes;
TypeID = find(strcmp(Types(:,4),'jpeg'));
setappdata(hObject, 'FileMenuFcnLastExportedAsType', TypeID)

% % Change the default filename (I presume it is not used otherwise) => causes troubles for GUIDE
% [Dum FigName] = fileparts(get(hObject, 'filename'));
% set(hObject, 'filename', fullfile(pwd, [FigName Types{TypeID,3}]))


% --- Outputs from this function are returned to the command line.
function varargout = dd_patch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isfield(handles, 'close')
    H = handles.axes_SumRes;
    clear handles                % Release working memory (?)
    waitfor(H);
    handles = guidata(hObject);  % Retrieve latest handles
end

% Return if there was nothing to process (i.e. no DTIFiles)
if ~isfield(handles, 'DWISel')
	varargout{1} = [];
	varargout{2} = [];
	if ishandle(hObject)
		close(hObject)
	end
	return
end

% Find out which are the selected outliers
DWISelI       = find(handles.DWISel);
[Regqi Regzi] = ind2sub(size(handles.Vol(:,:,1,1)), find(handles.Outliers==2));
[Slcqi Slczi] = ind2sub(size(handles.Vol(:,:,1,1)), find(handles.Outliers==3));
Output        = prepend(handles.input, 'artc_', DWISelI(unique([Regqi' Slcqi'])));

% Correct the voxel and slicewise artefacts
for n = 1:numel(Regqi)
    handles.Vol(Regqi(n), Regzi(n), handles.PCMask(Regqi(n), Regzi(n), :)) = ...
        handles.CVol(Regqi(n), Regzi(n), handles.PCMask(Regqi(n), Regzi(n), :));
end
for n = 1:numel(Slcqi)
    handles.Vol(Slcqi(n), Slczi(n), :) = handles.CVol(Slcqi(n), Slczi(n), :);
end

% Print and save the corrected data & converged voxel and slice weightings
HG = spm_figure('FindWin', 'Graphics');
if ishandle(HG)
	LogName  = getappdata(HG, 'LogName');
	SeriesNr = getappdata(HG, 'SeriesNr');
	FIDLog	 = fopen([LogName(1:end-2) 'txt'], 'a');
	if FIDLog > -1
		fprintf(FIDLog, 'S%g\tPATCH:\tVoxel outliers =\t%g\tSlice outliers =\t%g\tNr of slices =\t%g\n', ...
			SeriesNr, sum(handles.Outliers==2), sum(handles.Outliers==3), numel(handles.Outliers));
		fclose(FIDLog);
	else
		warning('DIDI:PATCH:OpenFile', 'Could not open: %s', [LogName(1:end-2) 'txt'])
	end
else
	SeriesNr = 1;
end
myspm_print(sprintf('S%g: PATCH - estimated weights', SeriesNr), fileparts(handles.input(1,:)))
fprintf('------------------------\n')
fprintf('Accepted voxel outliers:  %g/%g\n', sum(handles.Outliers==2), numel(handles.Outliers))
fprintf('Accepted slicewise outliers: %g/%g\n', sum(handles.Outliers==3), numel(handles.Outliers))
fprintf('Saving %g artefact-corrected files\n', numel(unique([Regqi' Slcqi'])))
for q = unique([Regqi' Slcqi'])
    if numel(q)==0, break, end   % q = 0-by-1 array
    Hdr       = spm_vol(handles.input(DWISelI(q),:));
	Hdr       = rmfield(Hdr, 'pinfo');
    Hdr.fname = Output(DWISelI(q),:);
    spm_write_vol(Hdr, ipermute(handles.Vol(q,:,:,:), handles.tzyx));
    [Pth Nme] = fileparts(handles.input(DWISelI(q),:));
	if exist(fullfile(Pth, ['artc_' Nme '.mat']), 'file')
		delete(fullfile(Pth, ['artc_' Nme '.mat']))	% copyfile can't overwrite writable files that it doesn't own
	end
    copyfile(fullfile(Pth, [Nme '.mat']), fullfile(Pth, ['artc_' Nme '.mat']))
end

% Save the WLS weights
PATCH.input    = handles.input;
PATCH.MaskFile = handles.MaskFile;
PATCH.b0File   = handles.b0File;
PATCH.WCrit    = handles.WCrit;
PATCH.WVox     = handles.WVox;
PATCH.WSlc     = handles.WSlc;
PATCH.WLog     = handles.WLog;
PATCH.Thresh   = handles.Thresh;
PATCH.Outliers = handles.Outliers;
PATCH.AffM     = handles.AffM;
disp(['Saving: ' fullfile(fileparts(handles.input(1,:)), 'PATCH')])
save(fullfile(fileparts(handles.input(1,:)), 'PATCH'), 'PATCH')

% Save the sum of the total weights as a nifti-file
Sz		  = size(handles.WVox);
Ext 	  = ['.' getfield(spm_get_defaults('images'),'format')];
WT		  = handles.Mask .* shiftdim(sum(1 - bsxfun(@times, handles.WVox, handles.WSlc)));
Hdr		  = spm_vol(handles.input(1,:));
Hdr       = rmfield(Hdr, 'pinfo');
Hdr.fname = fullfile(fileparts(Hdr.fname), ['PATCH_W' Ext]);
Hdr.dt(1) = spm_type('float32');
disp(['Saving: ' Hdr.fname])
spm_write_vol(Hdr, ipermute(WT, handles.tzyx(2:end)));  % NB the 4th dim is always t

% Save the nr of outliers as a nifti-file
WT = zeros(Sz);
for n = 1:numel(Regqi)
    WT(Regqi(n), Regzi(n), handles.PCMask(Regqi(n), Regzi(n), :)) = 1;
end
for n = 1:numel(Slcqi)
    WT(Slcqi(n), Slczi(n), :) = 1;
end
WT        = shiftdim(sum(WT)) .* handles.Mask;
Hdr       = spm_vol(handles.input(1,:));
Hdr       = rmfield(Hdr, 'pinfo');
Hdr.fname = fullfile(fileparts(Hdr.fname), ['PATCH_N' Ext]);
disp(['Saving: ' Hdr.fname])
spm_write_vol(Hdr, ipermute(WT, handles.tzyx(2:end)));  % NB the 4th dim is always t

% Finish up
varargout{1} = Output;
varargout{2} = PATCH;
if ishandle(hObject)
	close(hObject)
end


% --- Executes on key press with focus on figure and of its controls.
function dd_patch_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to dd_patch (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if length(handles.Outliers)<=1 && ~isfield(handles, 'close')
    return
end

switch [eventdata.Modifier{:} eventdata.Key]
    case {'1' 'delete'}
        radiobutton_DWI_Callback(handles.radiobutton_DWI, [], handles)
    case {'2' 'insert'}
        radiobutton_PCDWI_Callback(handles.radiobutton_PCDWI, [], handles)
    case '3'
        radiobutton_CDWI_Callback(handles.radiobutton_CDWI, [], handles)
    case {'leftarrow' 'backspace'}
        Ind = get(handles.slider_Ind, 'Value');
        if Ind > 1
            set(handles.slider_Ind, 'Value', Ind - 1)
            slider_Ind_Callback(handles.slider_Ind, [], handles)
        end
    case {'rightarrow' 'space'}
        Ind = get(handles.slider_Ind, 'Value');
        if Ind < get(handles.slider_Ind, 'Max')
            set(handles.slider_Ind, 'Value', Ind + 1)
            slider_Ind_Callback(handles.slider_Ind, [], handles)
        end
    case {'uparrow' 'controluparrow' 'altuparrow' 'shiftuparrow'}
        pushbutton_zplus_Callback(handles.pushbutton_zplus, [], handles)
    case {'downarrow' 'controldownarrow' 'altdownarrow' 'shiftdownarrow'}
        pushbutton_zmin_Callback(handles.pushbutton_zmin, [], handles)
    case {'controlrightarrow' 'altrightarrow' 'shiftrightarrow'}
        pushbutton_qplus_Callback(handles.pushbutton_qplus, [], handles)
    case {'controlleftarrow' 'altleftarrow' 'shiftleftarrow'}
        pushbutton_qmin_Callback(handles.pushbutton_qmin, [], handles)
    case 'escape'
        set(handles.slider_Ind, 'Value', get(handles.slider_Ind, 'Value'))
        slider_Ind_Callback(handles.slider_Ind, [], handles)
    case 'home'
        set(handles.slider_Ind, 'Value', 1)
        slider_Ind_Callback(handles.slider_Ind, [], handles)
    case 'end'
        set(handles.slider_Ind, 'Value', get(handles.slider_Ind, 'Max'))
        slider_Ind_Callback(handles.slider_Ind, [], handles)
    case {'controlreturn' 'controlenter'}
        pushbutton_accept_Callback(handles.pushbutton_accept, [], handles)
    case 'f1'
        doc('dd_patch')
end


% --- Executes on button pr1ess in radiobutton_DWI.
function radiobutton_DWI_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_DWI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Outliers(handles.Ind) = 1;
set(handles.radiobutton_DWI, 'Value', 1)
set(handles.radiobutton_PCDWI, 'Value', 0)
set(handles.radiobutton_CDWI, 'Value', 0)
guidata(hObject, handles);
set(findobj(handles.axes_DWI, 'Tag', 'DWIBorder'), 'Visible', 'on')
set(findobj(handles.axes_PCDWI, 'Tag', 'PCDWIBorder'), 'Visible', 'off')
set(findobj(handles.axes_CDWI, 'Tag', 'CDWIBorder'), 'Visible', 'off')

resimplot(handles, true)


% --- Executes on button press in radiobutton_PCDWI.
function radiobutton_PCDWI_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_PCDWI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setslider(handles, 2)

handles.Outliers(handles.Ind) = 2;
set(handles.radiobutton_DWI, 'Value', 0)
set(handles.radiobutton_PCDWI, 'Value', 1)
set(handles.radiobutton_CDWI, 'Value', 0)
guidata(hObject, handles);
set(findobj(handles.axes_DWI, 'Tag', 'DWIBorder'), 'Visible', 'off')
set(findobj(handles.axes_PCDWI, 'Tag', 'PCDWIBorder'), 'Visible', 'on')
set(findobj(handles.axes_CDWI, 'Tag', 'CDWIBorder'), 'Visible', 'off')

resimplot(handles, true)


% --- Executes on button press in radiobutton_CDWI.
function radiobutton_CDWI_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_CDWI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setslider(handles, 3)

handles.Outliers(handles.Ind) = 3;
set(handles.radiobutton_DWI, 'Value', 0)
set(handles.radiobutton_PCDWI, 'Value', 0)
set(handles.radiobutton_CDWI, 'Value', 1)
guidata(hObject, handles);
set(findobj(handles.axes_DWI, 'Tag', 'DWIBorder'), 'Visible', 'off')
set(findobj(handles.axes_PCDWI, 'Tag', 'PCDWIBorder'), 'Visible', 'off')
set(findobj(handles.axes_CDWI, 'Tag', 'CDWIBorder'), 'Visible', 'on')

resimplot(handles, true)


% --- Executes on slider movement.
function slider_Ind_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Ind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.slider_Ind, 'Visible'), 'off')
    initimpanels(handles, 'On')
end

Ind = round(get(hObject,'Value'));  % Position of slider
set(hObject, 'Value', Ind)          % Use round numbers only (> mouse interaction)

Outliers    = find(handles.Outliers);
handles.Ind = Outliers(Ind);
guidata(hObject, handles);

resimplot(handles)


% --- Executes during object creation, after setting all properties.
function slider_Ind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Ind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_Help.
function pushbutton_Help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

doc(mfilename)


% --- Executes on button press in pushbutton_Reset.
function pushbutton_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Outliers(:) = 0;   % Reset the PCDWI settings
handles.Outliers(handles.WMin(:) < handles.Thresh(1)) = 2;  % Set the new PCDWI settings
handles.Outliers(handles.WSlc(:) < handles.Thresh(2)) = 3;  % Set the new CDWI settings
fprintf('Redetected voxel outliers:  %g/%g\n', sum(handles.Outliers==2), numel(handles.WMin))
fprintf('Redetected slicewise outliers: %g/%g\n', sum(handles.Outliers==3), numel(handles.WMin))

% Re-initialize the panels
handles.Ind = 1;
NrOutliers  = numel(find(handles.Outliers));
if NrOutliers == 0
    initimpanels(handles, 'Off')    % Make image-panels invisible
end
if NrOutliers <= 1                  % There is outlier slider to show
    set(handles.slider_Ind, 'Value', 1, 'Min', 1-eps, 'Max', 1, 'SliderStep', [1 1])
else
    set(handles.slider_Ind, 'Value', 1, 'Min', 1, 'Max', NrOutliers, ...
        'SliderStep', [1 1]/(NrOutliers-1))
end
guidata(hObject, handles);

% Draw first slice
if NrOutliers >= 1
    slider_Ind_Callback(handles.slider_Ind, [], handles)
end


% --- Executes on button press in pushbutton_Threshold.
function pushbutton_Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the mouse input and set the threshold
[Dum y] = ginput(1);
YLim	= get(gca, 'YLim');
DWI		= handles.Outliers==1;					% Store the DWI settings
if strcmp(get(gca, 'Tag'), 'axes_MinRes') && y>=YLim(1) && y<=YLim(2)
    handles.Thresh(1) = y;
    CDWI = handles.Outliers==3;					% Store the CDWI settings
    handles.Outliers(handles.Outliers==2) = 0;	% Reset the PCDWI settings
    handles.Outliers(handles.WMin(:) < y) = 2;	% Set the new PCDWI settings
    handles.Outliers(CDWI)                = 3;	% Restore the CDWI settings
    set(findobj(handles.axes_MinRes, 'Tag', 'ThresLijn'), 'YData', [y y])
    fprintf('New min(slice) threshold: %g\n', y)
    handles.PCMask = handles.WVox < y;			% Update the PCDWI-mask
    handles.PCMask(:,~handles.Mask) = false;
elseif strcmp(get(gca, 'Tag'), 'axes_SumRes') && y>=YLim(1) && y<=YLim(2)
    handles.Thresh(2) = y;
    PCDWI = handles.Outliers==2;				% Store the PCDWI settings
    Combi = handles.Outliers==3 & handles.WMin(:)<handles.Thresh(1);
    handles.Outliers(handles.Outliers==3) = 0;	% Reset the CDWI settings
    handles.Outliers(PCDWI | Combi)       = 2;	% Restore the PCDWI settings
    handles.Outliers(handles.WSlc(:) < y) = 3;	% Set the new CDWI settings
    set(findobj(handles.axes_SumRes, 'Tag', 'ThresLijn'), 'YData', [y y])
    fprintf('New sum(slice) threshold: %g\n', y)
else
    return
end
handles.Outliers(DWI) = 1;						% Restore the DWI settings
fprintf('Redetected voxel outliers:  %g/%g\n', sum(handles.Outliers==2), numel(handles.WMin))
fprintf('Redetected slicewise outliers: %g/%g\n', sum(handles.Outliers==3), numel(handles.WMin))

% Re-initialize the panels
handles.Ind = 1;
NrOutliers  = numel(find(handles.Outliers));
if NrOutliers == 0
    initimpanels(handles, 'Off')    % Make image-panels invisible
end
if NrOutliers <= 1                  % There is outlier slider to show
    set(handles.slider_Ind, 'Value', 1, 'Min', 1-eps, 'Max', 1, 'SliderStep', [1 1])
else
    set(handles.slider_Ind, 'Value', 1, 'Min', 1, 'Max', NrOutliers, ...
        'SliderStep', [1 1]/(NrOutliers-1))
end
guidata(hObject, handles);

% Draw first slice
if NrOutliers >= 1
    slider_Ind_Callback(handles.slider_Ind, [], handles)
end


% --- Executes on button press in pushbutton_recalc.
function pushbutton_recalc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_recalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

WCrit = [];
while numel(WCrit)<1 || numel(WCrit)>5
    Ans = inputdlg({'Weighting criteria:' 'Re-use weighting parameters'}, ...
                   'IRLS recalculation', 1, {sprintf('%g  %g  %g  %g  %g',handles.WCrit) 'No'});
    if isempty(Ans)
        return
    else
        WCrit = str2num(Ans{1});
    end
end
handles.WCrit = WCrit;
if ~strcmp(Ans{2}, 'Yes')
    handles.WVox = [];
    handles.Slc  = [];
end
set(hObject, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
set(handles.pushbutton_accept, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
drawnow
dd_patch_OpeningFcn(handles.dd_patch, eventdata, handles, ...
    {handles.input}, handles.Mask, handles.WCrit, false, handles.AffM, handles.Outliers)
set(hObject, 'Enable', 'on', 'ForegroundColor', [0 0 0])
set(handles.pushbutton_accept, 'Enable', 'on', 'ForegroundColor', [0 0 0])


% --- Executes on button press in pushbutton_accept.
function pushbutton_accept_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stop user from hitting buttons
set(hObject, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
set(handles.pushbutton_recalc, 'Enable', 'inactive', 'ForegroundColor', [0.5 0.5 0.5])
drawnow

% Closing this figure resumes the outputfcn
delete(handles.axes_SumRes)


% --- Executes on button press in pushbutton_zplus.
function pushbutton_zplus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_zplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NewInd = handles.Ind + size(handles.Vol,1);
if NewInd <= length(handles.Outliers)
    handles.Ind = NewInd;
end
guidata(hObject, handles);

resimplot(handles)


% --- Executes on button press in pushbutton_zmin.
function pushbutton_zmin_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NewInd = handles.Ind - size(handles.Vol,1);
if NewInd >= 1
    handles.Ind = NewInd;
end
guidata(hObject, handles);

resimplot(handles)

    
% --- Executes on button press in pushbutton_qplus.
function pushbutton_qplus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_qplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~exist('handles', 'var')
    handles = guidata(hObject);
end
handles.Ind = min([handles.Ind+1 length(handles.Outliers)]);
guidata(hObject, handles);

resimplot(handles)


% --- Executes on button press in pushbutton_qmin.
function pushbutton_qmin_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_qmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Ind = max([handles.Ind-1 1]);
guidata(hObject, handles);

resimplot(handles)


% --- Executes on button press in pushbutton_Set.
function pushbutton_Set_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x Dum] = ginput(1);
XLim	= get(gca, 'XLim');
if (strcmp(get(gca, 'Tag'), 'axes_SumRes') || strcmp(get(gca, 'Tag'), 'axes_MinRes')) ...
        && x>=XLim(1) && x<=XLim(2)
    handles.Ind = round((x-1) * size(handles.Vol,1) + 1);
    guidata(hObject, handles);
    resimplot(handles)
end

%% ------------------------- END -----------------------


function [Vol WVox WSlc WLog SclF] = runpatch(Vol, VolHdr, Mask, b0Vol, q, WCrit, WVox, WSlc, WLog, Method, IsRecalc, AffM, FNames)

% This is the main function where all the calculations are done
%
% INPUT:
% Vol		- tzyx matrix of DWI-volumes (native space)
% VolHdr	- Header of the first DWI-volume
% Mask		- Logical zyx volume
% b0Vol		- tzyx matrix of bo-volume (reference space). Useful for catching IRLS failures
% q			- t*3 matrix with diffusion gradients
% WCRit		- a 5-element array with the weighting criteria
% WVox		- tzyx weighting volume from a previous run (can be empty)
% WSlc		- tz weighting volume from a previous run (can be empty)
% WLog		- tzyx weighting volume from a previous run (can be empty)
% Method	- String for method selection
% IsRecalc	- Boolean that is true for usage IV (= Recalculation from the GUI)
% AffM		- Realignment coordinate transformations in world coordinates (4x4xt)
% FNames	- DWI filenames (for displaying purpose)

% Initialize variables
MinIter  = 5;							% Minimum number of iterations
MaxIter  = 20;							% Maximum number of iterations
CThres	 = 0.01;						% Change threshold (=> fractional weights)
NThres	 = 500;							% Convergence threshold (=> tolerated nr of changes)
RThres	 = 50;							% Convergence threshold (=> tolerated nr of runaways)
%MedKrn  = double(getnhood(strel('disk', 3, 0)));	% Kernel for median filtering
MedKrn   = ones(5);						% Kernel for median filtering
%MorKrn  = strel('disk', 2, 0);			% Kernel for morphological operations
R		 = 4;							% kernel radius
MorKrn	 = ones(2*R+1, 2*R+1);			% Kernel for 2D morphological (dilation) operations
[I J]	 = ind2sub(size(MorKrn), find(MorKrn));
MorKrn(sqrt((I-R-1).^2 + (J-R-1).^2) > R) = 0;
R		 = 2;							% Radius of the 3D dilation/erosion kernel
Kernel   = ones(2*R+1, 2*R+1, 2*R+1);
[I J K]  = ind2sub(size(Kernel), find(Kernel));
Kernel(sqrt((I-R-1).^2 + (J-R-1).^2 + (K-R-1).^2) > R) = 0;
MaskD	 = logical(spm_dilate(double(Mask), Kernel));	% Use a larger mask in native space to catch motion (also needed when Mask is trimmed; not implemented, see below)
MaskE	 = logical(spm_erode(double(Mask), Kernel));	% Use a smalller mask for testing convergence
% Mask	 = logical(spm_erode(double(Mask), Kernel));	% Avoid rim effects, however, this does not work well when having large movements
[Zi Mi]  = find(Mask);					% Used in native space
ZiS      = unique(Zi)';
[ZiD MiD] = find(MaskD);				% Used in realigned space
% ZiS		 = unique(ZiD)';
VolSz    = size(Vol);
N0		 = numel(Zi) * VolSz(1);		% Nr of in-mask voxels * nr of volumes
Res		 = zeros(VolSz, class(Vol));	% Initialize Res volume
if ~isempty(AffM)
	Realign = true;
else
	Realign = false;
end
if nargin<7 || isempty(WVox)
	WVox  = ones(VolSz, class(Vol));	% Initialize weighting volume
	WSlc  = ones(VolSz(1:2));
	PrevW = false;
else
	PrevW = true;
end
WVoxSlc = ones(VolSz, class(Vol));		% Initialize WVoxSlc = WVox*WSlc
for n = 1:numel(Zi)
	WVoxSlc(:,Zi(n),Mi(n)) = WVox(:,Zi(n),Mi(n)) .* WSlc(:,Zi(n));
end
WVoxSlc(:,~Mask) = 1;
if strcmp(Method, 'PATCH') && WCrit(4)<1 && WCrit(5)<1
	VibrCorrect = true;
	Vibr		= abs(q(:,1)) > WCrit(5);	% The directions that may be corrupted
	fprintf('Nr of discarded x-gradient weighted images: %g/%g\n', sum(Vibr), numel(Vibr))
else
	VibrCorrect = false;
end
if WCrit(1)<0
	WCrit(1) = abs(WCrit(1));
	PosOutl	 = false;
else
	PosOutl  = true;
end

% Initialize the graphical output windows
HG = spm_figure('FindWin', 'Graphics');
HI = spm_figure('FindWin', 'Interactive');
if ishandle(HG)
	spm_figure('Clear', HG);
	spm('FigName', 'IRLS graphics', HG);
else
	HG = [];										% spm_figure returns an empty 0-by-1 matrix
end
if ishandle(HI)
	spm_figure('Clear', HI);
	spm('FigName', 'IRLS convergence', HI);
	HC  = axes('Position', [0.1 0.1 0.8 0.75], 'Parent',HI, 'Box','on');
	HW	= axes('Position', [0.1 0.85 0.8 0.1], 'Parent',HI);
	HWL = plot(HW, [0 0], [0.5 0.5], 'linewidth',10);
	HWT = text(0, 1, [Method ' 1:'], 'Parent', HW);
	set(HW, 'visible', 'off', 'xlim',[0 1], 'ylim',[0 1])
	set(HI, 'visible', 'on')
else
	HI=[]; HC=[]; HWL=[]; HWT=[];
end

% Compute the weight for the logarithmic distortion of the noise (= Basser et al., 1994)
if strcmp(Method, 'PATCH') && WCrit(3)>0
    LogCorrect = true;
	[Dum BG]   = dd_snr(Vol, Mask, FNames, 'Background noise (DWI)');	% Compute the background noise level
	BG0		   = BG;								% Store for later use
	if isempty(WLog)
		WLog = bsxfun(@rdivide, Vol, BG);			% The signal / mean background noise
		WLog = (1-WCrit(3)) + WCrit(3)*WLog;
		% NB: An easy and robust (but less sensitive) alternative would be to assume
		% that the signal noise is the same everywhere and just take WLog = Vol;
	end
else
	LogCorrect = false;
	WLog       = [];
end

% --- STRATEGY ---
% Create realigned volumes from which the residuals are computed. The residuals
% are then resliced to native space and used to estimate the weights, which are
% resliced back to realigned space for the next iteration.
% NB: Mask may not be exactly realigned with first DW volume!
if Realign
	for n = 1:size(q,1)
		% Rotate the gradients using the world-coordinate transformation
		q(n,:) = dd_rotategradients(q(n,:), AffM(:,:,n));

		% Convert AffM from world- to voxel-coordinates (see bireslice => spm_sample_vol)
		AffM(:,:,n) = VolHdr.mat \ AffM(:,:,n) * VolHdr.mat;
	end
	Vol0 = Vol;
	Vol  = bireslice(Vol, AffM, HWT);				% Volumes in realigned space
	if LogCorrect
		[Dum BG] = dd_snr(Vol);						% Compute the resliced background noise level
		BG		 = repmat(mean(BG(:)), size(BG));	% Assume a constant noise level (= more stable)
		WLog	 = Vol ./ mean(BG(:));				% Re-compute WLog (/BG is just for correctness but has no effect)
		WLog	 = (1-WCrit(3)) + WCrit(3)*WLog;
	end
end

% Calculate the residual forming matrix (t = [Txx Tyy Tzz Txy Txz Tyz]; y = Qt + e)
% [Dum DimFlip] = dd_rotategradients;
% q = repmat(DimFlip, size(q(:,1))) .* q;			% Flip to world coordinates is without consequense
Q = [q.^2 2*q(:,1).*q(:,2) 2*q(:,1).*q(:,3) 2*q(:,2).*q(:,3)];
I = eye(size(Q,1));

% Go to the log-domain and avoid log(0) = -Inf and log(-#) = complex
Mask(shiftdim(any(Vol<=0))) = false;				% Avoid voxels with outside volume sampling altogether
[Zi Mi]		= find(Mask);							% Used in native space
Vol(Vol<=0) = 1;
Vol			= log(Vol);
b0Vol		= log(b0Vol);

% Begin iterative re-weighted least squares regression (see also wlsqregress)
NrUpd        = round(numel(Zi)/250);				% Divide display waitbar in 250 pieces
MaxChange	 = [];
Changed		 = [];
Iterate      = true;
nIter        = 0;
CondNr		 = 10*eps(class(WVox));
disp('Starting IRLS estimation...')
while Iterate
    
    %--- Compute the weighted residuals (inside the mask only). NB: In realigned space
	nIter		= nIter + 1;
	Runaway		= 0;
	WVoxPrev	= WVox(:, MaskE);					% Remember the old weights (native sace)
	WVoxSlcPrev = WVoxSlc(:, MaskE);				% Remember the old weights (native sace)
	if nIter==1 && ~PrevW
		QiQ         = Q/(Q'*Q)*Q';					% Q times the pseudoinverse of Q
		Res(:,Mask) = (I - QiQ) * Vol(:,Mask);		% Simpler & faster when WT=ones
	else
		if Realign									% Reslice the weights again to realigned space
			WVoxSlc = exp(bireslice(log(WVoxSlc), AffM, HWT));	% Nonlinearity of W is undesired when reslicing
		end
        for n = 1:numel(Zi)							% Construct the voxelwise total-weights
			WT = WVoxSlc(:,Zi(n),Mi(n));			% The voxel + slicewise weights
            if LogCorrect
                WT = WT .* WLog(:,Zi(n),Mi(n));		% Add weight for log distortion of the noise
            end
			if VibrCorrect
				WT(Vibr) = WT(Vibr) * WCrit(4);		% Add weight for x-vibrations
			end
			WT(Vol(:,Zi(n),Mi(n))==0) = 0;			% Ignore missing values (i.e. zeros)
			WQ = bsxfun(@times, WT, Q);				% Faster than repmat
			if rcond(WQ'*WQ) > CondNr				% Make sure the matrix inverse can be properly computed
%			if sum(WT>0) >= 6						% Make sure we have enough data to make a new estimate
				QiWQW = Q/(WQ'*WQ)*WQ'*diag(WT);
				Res(:,Zi(n),Mi(n)) = (I - QiWQW) * Vol(:,Zi(n),Mi(n));
			end
			if any(abs(Res(:,Zi(n),Mi(n))) > b0Vol(Zi(n),Mi(n))) % Catch extreme residuals (i.e. obvious IRLS failures)
				[Yi Xi] = ind2sub(VolSz(3:4), Mi(n));
				fprintf('Resetting runaway IRLS regression [z=%d, y=%d, x=%d]\n', Zi(n), Yi, Xi)
				Res(:,Zi(n),Mi(n)) = 0;				% This will reset all weights to 1 and preserve the original data
				Runaway = Runaway + 1;
			end
            if rem(n, NrUpd)==0
				mywaitbar(1-n/numel(Zi), HWL, sprintf('%s %g: Computing residuals (n=%g)', Method, nIter, numel(Zi)-n), HWT)
            end
        end
	end
	Res0 = Res;										% Keep an unfiltered copy (see end of iteration)
	
	%--- Compute the new weights. NB: In native space
	mywaitbar(0, HWL, sprintf('%s %g: Computing weights (z=%g)', Method, nIter, ZiS(1)), HWT)
	switch Method
	 case 'PATCH'
		
		% Subtract the median residual for potentially faster and more robust convergence (= PATCH+)
		Res(:,Mask) = bsxfun(@minus, Res(:,Mask), median(Res(:,Mask)));
		
		% Compute the new logarithmic distortion weight using the estimated signal (= Salvador et al., 2005)
		% if LogCorrect && ~VibrCorrect						% Using the measured values helps to fight the vibation artefacts
		if LogCorrect && ~VibrCorrect						% One iteration should already be enough according to Salvador; this could possibly be more stable
			WLog = bsxfun(@rdivide, exp(Vol - Res), BG);	% exp(Vol-Res) = estimated signal
			WLog = (1-WCrit(3)) + WCrit(3)*WLog;			% NB: In realigned space
		end

		if Realign											% Reslice the residuals back to native space
			Res = bireslice(Res, AffM, HWT, 'inv');
		end
		
		% Compute the slicewise weights
		ZiSSel = ZiS(sum(MaskE(ZiS,:),2) > 4*sum(MorKrn(:)));	% Make sure we have more power than WVox
		for zi = ZiSSel
			WSlc(:,zi) = sum(Res(:, zi, MaskE(zi,:)), 3) / sqrt(sum(MaskE(zi,:)));
		end
		SclF = abs(WSlc(:,ZiSSel));
		SclF = 1.4826 * (median(SclF(:)) + median(SclF))/2;		% Scaling per slice + per volume
		WSlc(:,ZiSSel) = exp(-bsxfun(@rdivide, WCrit(2) * WSlc(:,ZiSSel), SclF).^2);

		% Compute the voxel weights
		SclF = 1.4826 * smoothn(shiftdim(median(abs(Res))), 2);
		if PosOutl                  % Process the residuals using median filtering (ordfilt2=medfilt2 with more flexible kernel)
			Resyxtz	= permute(Res, [3 4 1 2]);
			for zi = ZiS
			%     for qi = 1:VolSz(1)	% Old (slow) for loop
			% 		Res(qi,zi,:,:) = zmedfilt2s(shiftdim(Res(qi,zi,:,:)), single(MedKrn));
					% if rem(sum(MedKrn(:)),2) == 1	% The kernel has a odd number of ones
					% 	Res(qi,zi,:,:) = ordfilt2(shiftdim(Res(qi,zi,:,:)), (sum(MedKrn(:))+1)/2, MedKrn);
					% else
					% 	Res(qi,zi,:,:) = (ordfilt2(shiftdim(Res(qi,zi,:,:)), sum(MedKrn(:))/2, MedKrn) + ...
					% 					  ordfilt2(shiftdim(Res(qi,zi,:,:)), 1+sum(MedKrn(:))/2, MedKrn))/2;
					% end
			%     end
				Resyxtz(:,:,:,zi) = zmedfilt2sp(Resyxtz(:,:,:,zi), int8(MedKrn));
				mywaitbar(zi/ZiS(end), HWL, sprintf('%s %g: Computing weights (z=%g)', Method, nIter, zi), HWT)
			end
			Res = ipermute(Resyxtz, [3 4 1 2]);
		end
		WVox(:,MaskE) = exp(-bsxfun(@rdivide, WCrit(1) * Res(:,MaskE), SclF(MaskE)').^2);
		if ~PosOutl						% Consider negative outliers only for WVox (NB: this may bias the estimates in the absence of artefacts)
			WVox(Res>0)	= 1;
			WVoxyxtz	= permute(WVox, [3 4 1 2]);
			for zi = ZiS
			%     for qi = 1:VolSz(1)	% Old (slow) for-loop
			%         WVox(qi,zi,:,:) = spm_erode(spm_dilate(double(shiftdim(WVox(qi,zi,:,:))), MorKrn), MorKrn);
			%     end
				WVoxyxtz(:,:,:,zi) = zerode2sp(zdilate2sp(WVoxyxtz(:,:,:,zi), int8(MorKrn)), int8(MorKrn));
			end
			WVox = ipermute(WVoxyxtz, [3 4 1 2]);
		end
		% Expand (i.e. erode the weights) a little around the affected area
		WVoxyxtz = permute(WVox, [3 4 1 2]);
		for zi = ZiS
			% 	for qi = 1:VolSz(1)		% Old (slow) for-loop
			% 		% WVox(qi,zi,:,:) = imerode(shiftdim(WVox(qi,zi,:,:)), MorKrn);
			% 		WVox(qi,zi,:,:) = spm_erode(double(shiftdim(WVox(qi,zi,:,:))), MorKrn); % = ~10* slower than imerode but avoids use of toolbox
			% 	end
			WVoxyxtz(:,:,:,zi) = zerode2sp(WVoxyxtz(:,:,:,zi), int8(MorKrn));
			mywaitbar(1, HWL, sprintf('%s %g: Eroding weights (z=%g)', Method, nIter, zi), HWT)
		end
		WVox		   = ipermute(WVoxyxtz, [3 4 1 2]);
		WVox(:,~MaskE) = 1;

		% Combine the weights for the convergence evaluation and passage to the next iteration
		for n = 1:numel(Zi)				% NB: Zi indexes Mask, not MaskD
			WVoxSlc(:,Zi(n),Mi(n)) = WVox(:,Zi(n),Mi(n)) .* WSlc(:,Zi(n));
		end
		WVoxSlc(:,~MaskE) = 1;			% Multiple imerode operations can expand WVox far outside the mask
		
	 case 'GMM'
		 
		SclF = 1.4826 * shiftdim(median(abs(Res)));
		% WVox(:,Mask) = 1 ./ sqrt(SclF.^2 + Res(:,Mask).^2);							% = GMM
		WVox(:,Mask)	= 1 ./ sqrt(1 + bsxfun(@rdivide, Res(:,Mask), SclF(Mask)').^2);	% = GMM*SclF = scaled to 1 to allow easy thresholding (NB this does not effect the tensor estimate)
		WVoxSlc(:,Mask) = WVox(:,Mask);
		
	end
	
	%--- Determine if the weights have convergenced and display the covergence progres
	Changes			 = abs(WVoxSlc(:,MaskE) - WVoxSlcPrev) ./ WVoxSlcPrev;
	Changes(isinf(Changes)) = 0;
	MaxChange(end+1) = max(max(Changes));
	Changed(end+1)	 = sum(sum(Changes>CThres & WVoxSlc(:,MaskE)>CThres));	% Count the nr of voxels in which there was a significant change
	fprintf('Iteration %d:\tmax(|dW/W|) = %.3f\tChanged = %.4f (%d/%d)\n', nIter, MaxChange(end), Changed(end)/N0, Changed(end), N0)
	if ishandle(HC)
		plot(HC, [1:nIter]', [MaxChange(:) Changed(:)/N0 repmat(NThres/N0,[nIter 1])])
		legend(HC, {'max(|dW/W|)' 'Changed/N' 'Threshold'})
		xlabel(HC, 'Iteration nr')
		set(HC, 'XTick', unique(round(get(HC,'XTick'))), 'YScale','Log')				% 'YLim', [0 1])
	end
	plotweights(HG, WVox, WVoxPrev, WSlc, MaskE, nIter)
	Iterate = nIter<=MinIter | Changed(end)>NThres | Runaway>RThres;
	if (nIter >= MaxIter) || (PrevW && ~IsRecalc)	% PrevW && IsReCalc => Run the loop once when reload (Usage II)
		Iterate = false;
	end
	
	% Only update slices containing significant changes, i.e. make a new and smaller mask
	% In theory this is a great idea, but is a pain to implement as the new mask needs to
	% be resliced to realigned space and then dilated. It's probably not worth it.
	% NB: The code below is/was still WIP and left here as documentation
	% if Iterate && nIter > MinIter
	% 	[ZiE MiE] = find(MaskE);
	% 	ZiS		  = unique(ZiE)';
	% 	for zi = ZiS
	% 		if ~any(any(Changes(:,ZiE==zi) > CThres))
	% 			fprintf('Slice %d converged\n', zi)
	% 			Mask(zi,:)	= 0;
	% 			MaskD(zi,:) = 0;
	% 			MaskE(zi,:) = 0;
	% 		end
	% 	end
	% 	[Zi Mi]	  = find(Mask);
	% 	[ZiD MiD] = find(MaskD);
	% 	ZiS		  = unique(ZiD)';
	% end
	
end
fprintf('IRLS convergence:\tmax(|dW/W|) = %.3g\tChanged = %d\n', MaxChange(end), Changed(end))

%--- Compute the output variables [Vol WVox WSlc WLog (SclF)] in native space

% Check for runaway and non-converging regressions in the final results and give feedack to the user
if Changed(end)>NThres || Runaway>NThres			% Use NThres
	warning('DIDI:PATCH:Convergence', 'There were %d runaway and %d non-converging IRLS regressions in the data -- check your results if these numbers are very high', ...
			Runaway, Changed(end))
end

% Compute the predicted values in native space
Vol = exp(Vol - Res0);								% Realigned space
mywaitbar(1, HWL, sprintf('%s %g: Computing final corrected volumes', Method, nIter), HWT)
if Realign
	% Reslice the volumes back to native space (NB: WVox & WSlc still are)
	Vol = bireslice(Vol, AffM, HWT, 'inv');
	Vol(Vol==0) = Vol0(Vol==0);						% Account for out-of-volume sampling
end

% Compute WLog in native space
if LogCorrect
	WLog = bsxfun(@rdivide, Vol, BG0);
	WLog = (1-WCrit(3)) + WCrit(3) * WLog;			% Apply weighting criterium
end

if ishandle(HI)
	set(HI, 'visible', 'off')
end


function mywaitbar(x, HWL, String, HWT)

if ishandle(HWL)
	set(HWL, 'XData', [0 x])
	set(HWT, 'String', String)
	drawnow
end


function Vol = bireslice(Vol, AffMi, HWT, Inv)

if nargin<4, Inv=''; end
Inv	= strcmp(Inv, 'inv');
Txt = get(HWT, 'string');
Txt = [Txt(1:strfind(Txt, ':')) ' Reslicing volumes'];

VolSz	= size(Vol);						% Vol = tzyx
[Z Y X] = ndgrid(1:VolSz(2), 1:VolSz(3), 1:VolSz(4));	% Make z the fastest varying dimension (-> spm_sample_vol)
XYZ		= [X(:) Y(:) Z(:) ones(numel(X),1)]';
for ti = 1:VolSz(1)
	set(HWT, 'String', sprintf('%s (%g/%g)', Txt, ti, VolSz(1))), drawnow
	if Inv									% Inv = true => reslicing back to native space
		tXYZ = AffMi(:,:,ti) * XYZ;
	else									% Inv = false => reslicing to realigned space
		tXYZ = AffMi(:,:,ti) \ XYZ;			% NB: use inverse rotation in spm_sample_vol
	end
	% Reslice the xyz-volume using nn interpolation (trilinear interpolation didn't seem to work as well and is slower)
	Vol(ti,:) = spm_sample_vol(shiftdim(Vol(ti,:,:,:)), tXYZ(3,:), tXYZ(2,:), tXYZ(1,:), 0); % 1);
end


function plotweights(HG, WVox, WVoxPrev, WSlc, Mask, nIter)

if ~ishandle(HG)
	return
end
spm_figure('Clear', HG);
set(HG, 'visible', 'on')

mWVox	= mean(WVox);							% Average over volumes (time)
MinWVox	= min(mWVox(:));
if ~isfinite(MinWVox) || MinWVox==1
	MinWVox = 0;
end
h = subplot(5,2,1, 'Parent',HG);
imagesc(flipud(squeeze(min(mWVox,[],4))), 'Parent',h, [MinWVox 1])
axis(h, 'image')
ylabel(h, 'Voxel weights (MIP)')
title(h, sprintf('Iteration %g', nIter-1))
h = subplot(5,2,3, 'Parent',HG);
imagesc(flipud(squeeze(min(mWVox,[],3))), 'Parent',h, [MinWVox 1])
axis(h, 'image')
ylabel(h, 'Voxel weights (MIP)')
h = subplot(5,2,5, 'Parent',HG);
imagesc(flipud(squeeze(min(mWVox,[],2))), 'Parent',h, [MinWVox 1])
colorbar('peer', h)
axis(h, 'image')
ylabel(h, 'Voxel weights (MIP)')

dmWVox = zeros(size(mWVox));
dmWVox(1,Mask) = abs(mWVox(1,Mask) - mean(WVoxPrev));
MaxdWVox = max([dmWVox(:); eps]);
if ~isfinite(MaxdWVox)
	MaxdWVox = 1;
end
h = subplot(5,2,2, 'Parent',HG);	% Use an inverted grayscale
imagesc(-flipud(squeeze(max(dmWVox,[],4))), 'Parent',h, [-MaxdWVox 0])
axis(h, 'image')
title(h, sprintf('|Iteration %g - Iteration %g|', nIter-1, nIter-2))
h = subplot(5,2,4, 'Parent',HG);
imagesc(-flipud(squeeze(max(dmWVox,[],3))), 'Parent',h, [-MaxdWVox 0])
axis(h, 'image')
h = subplot(5,2,6, 'Parent',HG);
imagesc(-flipud(squeeze(max(dmWVox,[],2))), 'Parent',h, [-MaxdWVox 0])
HCB = colorbar('peer', h);
set(HCB, 'YDir','reverse', 'YTickLabel',num2str(abs(get(HCB,'YTick')')))
axis(h, 'image')

h = subplot(5,2,[7 8], 'Parent',HG);
plot(h, min(WVox(:,:,:), [], 3)')
hold(h, 'on')
plot(h, get(h, 'XLim'), 0.5*[1;1], 'k--')
axis(h, [1 size(WSlc,2) 0 1])
ylabel(h, 'Voxel weights')

h = subplot(5,2,[9 10], 'Parent',HG);
plot(h, WSlc')
hold(h, 'on')
plot(h, get(h, 'XLim'), 0.5*[1;1], 'k--')
axis(h, [1 size(WSlc,2) 0 1])
xlabel(h, 'Slice number (z)')
ylabel(h, 'Slice weights')


function wplot(hW, W, Thresh, Color, Tag)

% Plot the residual error distribution
[Nq Nz] = size(W);
hold(hW, 'on'), cla(hW)
plot(hW, 1:1/Nq:Nz+1-1/Nq, W(:), Color)
xlim(hW, [1 Nz+1])
ylim(hW, [0 1])
vertlijn(hW, 1,':', 'Tag', 'ZLijn')
plot(hW, 1,W(1),'xk', 1,W(1),'+y', 'Tag', 'ZMarker')
horzlijn(hW, Thresh, ['--' Color], 'Tag', 'ThresLijn')
xlabel(hW, 'Slice number (z)')
text((Nz+1)*0.98, 0.12, Tag, 'EdgeColor','k', 'HorizontalAlignment','right', 'Parent',hW)


function resimplot(handles, Hist)

%Kernel = strel('disk', 1);		% Used for making the ROI-outline
%Kernel = strel('square', 3);	% Used for making the ROI-outline
Kernel = ones(3);				% Used for making the ROI-outline

% Plot the selected slice
NqNz     = size(handles.Vol(:,:,1,1));
[qi zi]  = ind2sub(NqNz, handles.Ind);
ResM     = find(handles.Mask(zi, :, :));
WVox     = shiftdim(handles.WVox(qi, zi, :, :));
SclF     = shiftdim(handles.SclF(zi, :, :));
PCMask   = shiftdim(handles.PCMask(qi, zi, :, :));
DWI      = shiftdim(handles.Vol(qi, zi, :, :));
CDWI     = shiftdim(handles.CVol(qi, zi, :, :));
Res      = DWI - CDWI;
PCDWI    = DWI;
PCDWI(PCMask) = CDWI(PCMask);
DWIRange = [min([DWI(:); CDWI(:)]) max([DWI(:); CDWI(CDWI<5000)])];	% max(CDWI) can be very large if misestimated
YRange   = handles.BBox(2):-1:handles.BBox(1);						% Flip up-down
XRange   = handles.BBox(3):handles.BBox(4);
BigMask       = spm_dilate(double(PCMask(YRange, XRange)), Kernel);	% imdilate is faster but avoid use of toolbox
[EdgeX EdgeY] = findline(xor(BigMask, spm_erode(BigMask, Kernel)));
%[EdgeX EdgeY] = findline(edge(uint8(spm_dilate(PCMask(YRange, XRange), Kernel))));
%[EdgeX EdgeY] = findline(edge(uint8(PCMask(YRange, XRange))));
Border  = 1 + [0 XRange(end)-XRange(1) XRange(end)-XRange(1) 0 0
               0 0 YRange(1)-YRange(end) YRange(1)-YRange(end) 0];
ResRange = max(abs(Res(:)))*[-1 1];
if ResRange(1)==0
    ResRange = [-eps eps];
end
BVec = handles.q(qi,:);
if ~isempty(handles.AffM)
	[BVec DimFlip] = dd_rotategradients(BVec, inv(handles.AffM(:,:,qi)));	% Rotate BVec to voxel coordinates
else
	[Dum DimFlip] = dd_rotategradients;			% Just get the DimFlip
end
% Account for Volkmar's reference frame error; this should one day be replaced by a global default setting
BVec = DimFlip .* BVec;

DWIFiles	   = handles.input(handles.DWISel, :);
[Dum FileName] = fileparts(DWIFiles(qi,:));
set(handles.text_File, 'String', {FileName})
set(handles.text_Info, 'String', {sprintf('BVec = [%.2f %.2f %.2f];  SliceNr = %d', BVec, zi)})

%-- Plot axes
if handles.Outliers(handles.Ind) == 3
    hist(handles.axes_Hist, Res(ResM), 50)
    vertlijn(handles.axes_Hist, 0, '-k', 'Parent', handles.axes_Hist)
	vertlijn(handles.axes_Hist, mean(Res(ResM)), ':c')
    xlabel(handles.axes_Hist, 'Residual error')
    ylabel(handles.axes_Hist, 'Nr of voxels')
else
    imagesc(WVox(YRange, XRange), 'Parent', handles.axes_Hist, [0 1])
    set(handles.axes_Hist, 'XTick', [], 'YTick', [])
	axis(handles.axes_Hist, 'image')
    xlabel(handles.axes_Hist, 'Voxel weights')
end
if nargin>1 && Hist
    return
end

quiver(handles.axes_BVec, 0, 0, BVec(1), BVec(2), 'AutoScale', 'Off')
hold(handles.axes_BVec, 'on')
t = 0:2*pi/180:2*pi;
plot(handles.axes_BVec, sin(t), cos(t))
hold(handles.axes_BVec, 'off')
set(handles.axes_BVec, 'XLim', [-1 1], 'YLim', [-1 1], 'Visible', 'Off')

H_DWI = imagesc(DWI(YRange, XRange), 'Parent', handles.axes_DWI, DWIRange);
hold(handles.axes_DWI, 'on')
for n = 1:numel(EdgeX)
    plot(handles.axes_DWI, EdgeY{n}, EdgeX{n}, '-c', 'HitTest', 'Off')
end
plot(handles.axes_DWI, Border(1,:), Border(2,:), 'g', 'Tag','DWIBorder', 'Visible','Off')
hold(handles.axes_DWI, 'off')
axis(handles.axes_DWI, 'image')
colormap(handles.axes_DWI, 'gray')
colorbar('peer', handles.axes_DWI, 'Location', 'WestOutside')
set(handles.axes_DWI, 'XTick', [], 'YTick', [])
set(H_DWI, 'ButtonDownFcn', @(hObject, eventdata)radiobutton_DWI_Callback(H_DWI,[],handles))

H_PCDWI = imagesc(PCDWI(YRange, XRange), 'Parent', handles.axes_PCDWI, DWIRange);
hold(handles.axes_PCDWI, 'on')
plot(handles.axes_PCDWI, Border(1,:), Border(2,:), 'g', 'Tag','PCDWIBorder', 'Visible','Off')
hold(handles.axes_PCDWI, 'off')
axis(handles.axes_PCDWI, 'image')
set(handles.axes_PCDWI, 'XTick', [], 'YTick', [])
set(H_PCDWI, 'ButtonDownFcn', @(hObject, eventdata)radiobutton_PCDWI_Callback(H_PCDWI,[],handles))

H_CDWI = imagesc(CDWI(YRange, XRange), 'Parent', handles.axes_CDWI, DWIRange);
hold(handles.axes_CDWI, 'on')
plot(handles.axes_CDWI, Border(1,:), Border(2,:), 'g', 'Tag','CDWIBorder', 'Visible','Off')
hold(handles.axes_CDWI, 'off')
axis(handles.axes_CDWI, 'image')
set(handles.axes_CDWI, 'XTick', [], 'YTick', [])
set(H_CDWI, 'ButtonDownFcn', @(hObject, eventdata)radiobutton_CDWI_Callback(H_CDWI,[],handles))

imagesc(Res(YRange, XRange), 'Parent', handles.axes_Res, ResRange)
hold(handles.axes_Res, 'on')
for n = 1:numel(EdgeX)
    plot(handles.axes_Res, EdgeY{n}, EdgeX{n}, '-c')
end
hold(handles.axes_Res, 'off')
axis(handles.axes_Res, 'image')
colorbar('peer', handles.axes_Res, 'Location', 'WestOutside')
set(handles.axes_Res, 'XTick', [], 'YTick', [])
xlabel(handles.axes_Res, 'Residual errors')

imagesc(SclF(YRange, XRange), 'Parent', handles.axes_PC)
hold(handles.axes_PC, 'on')
for n = 1:numel(EdgeX)
    plot(handles.axes_PC, EdgeY{n}, EdgeX{n}, '-c')
end
hold(handles.axes_PC, 'off')
axis(handles.axes_PC, 'image')
set(handles.axes_PC, 'XTick', [], 'YTick', [])
xlabel(handles.axes_PC, 'Error variance')

% Move the vertlijn in the Res-figure
Count = 1;
while ~ishandle(handles.axes_SumRes) % Might happen when pc is slow and user clicks real fast?
    disp('Waiting for figure update...')
    drawnow  % This does not seem to help (but should)
    pause(0.01)
    Count = Count + 1;
    if Count > 250
        warning('DIDI:PATCH:GUIUpdate', 'Figure update is too slow')
        break
    end
end
HZ = findobj(handles.axes_MinRes, 'Tag', 'ZLijn');
set(HZ, 'XData', (zi + (qi-1)/NqNz(1))*[1 1])
HZ = findobj(handles.axes_MinRes, 'Tag', 'ZMarker');
set(HZ(:), 'XData', zi + (qi-1)/NqNz(1), 'YData', handles.WMin(qi,zi))
HZ = findobj(handles.axes_SumRes, 'Tag', 'ZLijn');
set(HZ, 'XData', (zi + (qi-1)/NqNz(1))*[1 1])
HZ = findobj(handles.axes_SumRes, 'Tag', 'ZMarker');
set(HZ(:), 'XData', zi + (qi-1)/NqNz(1), 'YData', handles.WSlc(qi,zi))

% Make sure the radio-buttons are visible
if strcmp(get(handles.radiobutton_DWI, 'Visible'), 'off')
    set(handles.radiobutton_DWI, 'Visible', 'On')
    set(handles.radiobutton_PCDWI, 'Visible', 'On')
    set(handles.radiobutton_CDWI, 'Visible', 'On')
end

% Set the DWI radiobuttons
if handles.Outliers(handles.Ind) <= 1
    set(handles.radiobutton_DWI, 'Value', 1)
    set(handles.radiobutton_PCDWI, 'Value', 0)
    set(handles.radiobutton_CDWI, 'Value', 0)
    set(findobj(handles.axes_DWI, 'Tag', 'DWIBorder'), 'Visible', 'on')
elseif handles.Outliers(handles.Ind) == 2
    set(handles.radiobutton_DWI, 'Value', 0)
    set(handles.radiobutton_PCDWI, 'Value', 1)
    set(handles.radiobutton_CDWI, 'Value', 0)
    set(findobj(handles.axes_PCDWI, 'Tag', 'PCDWIBorder'), 'Visible', 'on')
elseif handles.Outliers(handles.Ind) == 3
    set(handles.radiobutton_DWI, 'Value', 0)
    set(handles.radiobutton_PCDWI, 'Value', 0)
    set(handles.radiobutton_CDWI, 'Value', 1)
    set(findobj(handles.axes_CDWI, 'Tag', 'CDWIBorder'), 'Visible', 'on')
end

%---- Display WIP figure ----
try    wip = evalin('base', 'wip');
catch, wip = false; end
if wip && exist('dd_patch_wip.m', 'file')
    persistent PhaseVol
    dd_patch_wip
end


function setslider(handles, OutlierI)

% Expands the slider when undetected outliers are added
if ~handles.Outliers(handles.Ind)
    handles.Outliers(handles.Ind) = OutlierI;
    Outliers = find(handles.Outliers);
    IndVal   = find(Outliers==handles.Ind);
    if numel(Outliers)==1   % A rare situation (no outliers => set => use CDWI)
        Step = [1 1];
        set(handles.slider_Ind, 'Visible', 'On')
    else
        Step = [1 1]/(numel(Outliers)-1);
    end
    set(handles.slider_Ind, 'Value', IndVal, 'Max', numel(Outliers), ...
        'SliderStep', Step)
end


function initimpanels(handles, Tag)

set(handles.axes_DWI, 'Visible', 'Off')
set(handles.axes_PCDWI, 'Visible', 'Off')
set(handles.axes_CDWI, 'Visible', 'Off')
set(handles.axes_Res, 'Visible', 'Off')
set(handles.axes_PC, 'Visible', 'Off')
set(handles.axes_Hist, 'Visible', 'Off')

set(handles.radiobutton_DWI, 'Visible', Tag)
set(handles.radiobutton_PCDWI, 'Visible', Tag)
set(handles.radiobutton_CDWI, 'Visible', Tag)
set(handles.slider_Ind, 'Visible', Tag)
set(handles.OutlNav_text, 'Visible', Tag)

cla(handles.axes_DWI), colorbar('Off', 'peer', handles.axes_DWI)
cla(handles.axes_PCDWI)
cla(handles.axes_CDWI)
cla(handles.axes_Res), colorbar('Off', 'peer', handles.axes_Res)
cla(handles.axes_PC)
cla(handles.axes_Hist)
cla(handles.axes_BVec)


function PO = prepend(PI, pre, Sel)
%
% Input: PI  - filelist (SPM-style)
%        Sel - indexnr (not a logical)

if nargin<3
    Sel = 1:size(PI,1);
end
if isempty(Sel)
    Sel = false;
end

for n = 1:size(PI,1)
    [pth,nm,xt] = fileparts(PI(n,:));
    xt = strrep(xt, ',1', '');
    if any(Sel==n)
        PO{n} = fullfile(pth, [pre nm xt]);
    else
        PO{n} = fullfile(pth, [nm xt]);
    end
end
PO = char(PO);


function vertlijn(h, x, varargin)

% FUNCTION vertlijn(h, x, varargin)
%
% Plot een verticale lijn in je figuur op plaats 'x'

if nargin <= 2, varargin={'-'}; end
if isempty(x), return, end
x = x(:)';

n_as = get(h,'Nextplot');
hold(h, 'on')
y_lim = get(h,'YLim');
plot(h, [x;x], y_lim, varargin{:})
set(h,'Nextplot',n_as)


function horzlijn(h, y, varargin)

% FUNCTION horzlijn(h, y, <properties>)
%
% Plot een horizontale lijn in je figuur op hoogte 'y'

if nargin <= 2, varargin={'-'}; end
if isempty(y), return, end
y = y(:)';

n_as = get(h, 'Nextplot');
hold(h, 'on')
x_lim = get(h,'XLim');
plot(h, x_lim, [y;y], varargin{:})
set(h, 'Nextplot', n_as);


function [Xs Ys] = findline(XY)

% Orders X and Y such that they form a connected path

Xs = {[]};
Ys = {[]};

%[XY N] = bwlabel(XY, 4);							% ~10-50 times faster but uses a license :-(
[XY N] = roilabel2(XY, 4);

for n = 1:N

    [X Y] = find(XY==n);

    Xs{n} = X(1);
    Ys{n} = Y(1);
    while true
        [Dum I] = sort((X-X(1)).^2 + (Y-Y(1)).^2);	% Find the distances around the origin
        X		= X(I);								% Sort all points accordingly
        X(1)	= [];								% Remove the origin (which makes the nearest point the next origin)
        if isempty(X)
            Xs{n}(end+1) = Xs{n}(1);				% Close the path
            Ys{n}(end+1) = Ys{n}(1);
            break
        end
        Xs{n}(end+1) = X(1);						% Add the nearest point
        Y            = Y(I);
        Y(1)         = [];
        Ys{n}(end+1) = Y(1);
    end

end


function myspm_print(HdrTxt, FtrTxt)

% Robust against closed figures
HG = spm_figure('FindWin', 'Graphics');
if ishandle(HG)
	HD = axes('Position', [0.05 0.95 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0.75, 0.5, datestr(now), 'Parent',HD)
	text(0, 0.5, HdrTxt, 'Parent',HD)
	HD = axes('Position', [0.05 0 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0, 0.5, FtrTxt, 'Parent',HD)
	[S LWarn] = mywarning('Off','spm:spm_jobman:NotInitialised');
	spm_print(char(getappdata(HG, 'LogName')))
	mywarning(S, LWarn)
end


% ------------- copied from filemenufcn -----------------------------------
function list = localexporttypes

% build the list dynamically from printtables.m
[a,opt,ext,d,e,output,name] = printtables;                              %#ok

% only use those marked as export types (rather than print types)
% and also have a descriptive name
valid = strcmp(output,'X') & ~strcmp(name,'') & ~strcmp(d, 'QT'); 
name  = name(valid);
ext   = ext(valid);
opt   = opt(valid);

% remove eps formats except for the first one
iseps = strncmp(name,'EPS',3);
inds = find(iseps);
name(inds(2:end),:) = [];
ext(inds(2:end),:) = [];
opt(inds(2:end),:) = [];

for i=1:length(ext)
    ext{i} = ['.' ext{i}];
end
star_ext = ext;
for i=1:length(ext)
    star_ext{i} = ['*' ext{i}];
end
description = name;
for i=1:length(name)
    description{i} = [name{i} ' (*' ext{i} ')'];
end

% add fig file support to front of list
star_ext = {'*.fig',star_ext{:}};
description = {'MATLAB Figure (*.fig)',description{:}};
ext = {'.fig',ext{:}};
opt = {'fig',opt{:}};

[description,sortind] = sort(description);
star_ext = star_ext(sortind);
ext = ext(sortind);
opt = opt(sortind);

list = [star_ext(:), description(:), ext(:), opt(:)];
