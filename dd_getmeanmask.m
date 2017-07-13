function [Mask MeanImg Brain New] = dd_getmeanmask(Imgs, Realign, BETOpts, FName, SeriesNr)

% FUNCTION [Mask MeanImg Brain New] = DD_GETMEANMASK(Imgs, Realign, BETOpts, FName, SeriesNr)
%
% Get a previously made mask or otherwise make a new one with FSL BET.
%
% INPUT:
%	Imgs	 - filelist (SPM-style)
%	Realign  - If true then realign the images before taking the mean
%	BETOpts  - String array(s) that is passed to BET
%	FName	 - Name of the graphical output file. If FName is empty or NaN, then
%			   just get the requested files if they already exist and do nothing
%			   else, i.e. do not realign and do not save the displayed results.
%	SeriesNr - Series identifier that is only used in the graphical output
%
% OUTPUT:
%	Mask	 - String array
%
% Marcel, 24-4-2014
%
% See also: DD_BASICPROC_GETMASK

% Defaults
if nargin<5 || isempty(SeriesNr)
	SeriesNr = 1;
end
if nargin<4 || isempty(FName)
	FName = NaN;
end
if nargin<3 || isempty(BETOpts)
	BETOpts = '-R -m -f 0.35';
end
BETOpts = strtrim(BETOpts);
if nargin<2 || isempty(Realign)
	Realign = true;
end

% Get a graphics window
[S LWarn] = mywarning('Off', 'SPM:noDisplay');
spm_figure('GetWin', 'Graphics');
mywarning(S, LWarn)

% Filenames
Img1	= deblank(Imgs(1,:));
MeanImg = spm_file(Img1,	'prefix','mean');
Brain	= spm_file(MeanImg, 'suffix','_brain');
Mask	= spm_file(MeanImg, 'suffix','_brain_mask');

% Create or get the previously created 'mean'-image (= BET input image)
if ~exist(MeanImg,'file')
	if Realign && size(Imgs,1)>1
		% Use the first 10 images to create a mean image (just to speed things up a bit)
		Vols = spm_realign(Imgs(1:min(size(Imgs,1),10),:), struct('sep',2, 'rtm',1));
	else
		Vols = spm_vol(Imgs(1:min(size(Imgs,1),10),:));
	end
	spm_reslice(Vols, struct('mean',1, 'which',0));
end

% Get a brain mask
if exist(Mask,'file') && exist(Brain,'file')	% Reset orientations
	
	disp(['Using existing mask: ' Mask])
	spm_get_space(Mask,	   spm_get_space(MeanImg));
	spm_get_space(Brain,   spm_get_space(MeanImg));
	New		= false;
	BETOpts = '[PREVIOUS]';
	
else											% Create a new mask using BET	
	try
		
		Ext = spm_file(Img1,'ext');
		switch lower(Ext)
			case 'nii'
				Format = 'NIFTI';
			case 'img'
				Format = 'NIFTI_PAIR';
			otherwise
				error('Unknown output format .%s (%s)', Ext, Img1);
		end
		fprintf('\nCreating mask-file: %s\n', Mask);
		% THIS SHOULD WORK BUT BET DOES NOT RESPECT FSLOUTPUTTYPE [Sts Msg] = unix(['source ~/.bashrc; export FSLOUTPUTTYPE=' Format '; bet ' MeanImg ' ' Brain(1:end-4) ' ' BETOpts]);
		[Sts Msg] = unix(['source ~/.bashrc; bet ' MeanImg ' ' Brain(1:end-4) ' ' BETOpts]);
		if Sts % || ~isempty(Msg)
			if strfind(BETOpts, '-R')
				warning('DIDI:BET', 'Failure detected: Trying BET without the -R option')
				[Sts Msg] = unix(['source ~/.bashrc; bet ' MeanImg ' ' Brain(1:end-4) ' ' strrep(BETOpts,'-R','')]);
				if Sts % || ~isempty(Msg)
					error(Msg)
				end
			else
				error(Msg)
			end
		end
		unix(['source ~/.bashrc; fslchfiletype ' Format ' ' Mask ...
							  '; fslchfiletype ' Format ' ' Brain]);
		New	= true;
		
	catch BetError
		
		disp(['bet ' MeanImg ' ' Brain(1:end-4) ' ' BETOpts])
		warning('DIDI:BET', 'Automated mask creation on %s failed', MeanImg)
		disp(Msg)
		rethrow(BetError)
		
	end
end

% Display the results
iMask		= spm_vol(Mask);
iMask.fname = spm_file(tempname, 'ext','nii');
MaskVol		= spm_read_vols(spm_vol(Mask));
MeanVol		= spm_read_vols(spm_vol(MeanImg));
iMask       = spm_write_vol(iMask, ~MaskVol);				% Create a temporary inverted Mask
spm_check_registration(char(Brain, MeanImg))
spm_orthviews('Interp',0)									% NN interpolation is more crisp
spm_orthviews('context_menu','orientation',2)				% Use voxel-space of 1st image
spm_orthviews('AddColouredImage', 1, iMask, [0 1 1])		% [0 1 1] = 'cyan'
spm_orthviews('AddColouredImage', 2,  Mask,	[0 1 1])		% [1 0 0] = 'red'
spm_orthviews('XHairs','Off')
spm_orthviews('Window', 1:2, [0 0.5] * max(MeanVol(:)))		% Make the image a bit brighter
myspm_print(FName, sprintf('S%d: Brain-mask (bet %s)', SeriesNr, BETOpts))
delete(iMask.fname)


%% ------------------------- END ------------------------


function myspm_print(FileName, HdrTxt)

if isnan(FileName), return, end

% Robust against closed figures
HG = spm_figure('FindWin', 'Graphics');
if ishandle(HG)
	HL = findobj(HG, 'Tag','legend');
	for n = 1:numel(HL)
        set(HL(n), 'Box','Off')			% Legends scale badly on headless nodes
	end
	HD = axes('Position', [0.05 0.95 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0.75, 0.5, datestr(now), 'Parent',HD)
	text(0, 0.5, HdrTxt, 'Parent',HD)
	[S LWarn] = mywarning('Off','spm:spm_jobman:NotInitialised');
	spm_print(FileName)
	mywarning(S, LWarn)
end