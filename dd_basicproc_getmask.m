function [Mask Target] = dd_basicproc_getmask(Job, SubjNr, SeriesNr, LogName)

% FUNCTION [Mask Target] = DD_BASICPROC_GETMASK(Job, SubjNr, SeriesNr, LogName)
%
% Get a previously made mask or otherwise make a new one with BET. The DW images
% are realigned to the meanb0 image (if this Job option is set) and used to
% create the mask. Using the DW images works better than using the b0 images due
% to the suppressed CSF signal (so it cannot be mistaken for brain tissue).
%
% DD_BASICPROC_GETMASK is an internal function of dd_basicproc that is
% made available externally to allow distributed computing.
%
% INPUT:
%	SeriesNr - Series number(s) for which output images are returned (empty = all)
%	LogName  - If LogName is empty or NaN, then do nothing and just get the files,
% 			   i.e. do not realign the images and do not save the displayed results
%
% OUTPUT
%	Mask 	 - (cell)string array
%	Target	 - mean(bo-image)
%
% Marcel, 21-05-2014
%
% See also: DD_BASICPROC, DD_GETMEANMASK

% Defaults
if nargin<3 || isempty(SeriesNr)
	SeriesNr = 1:size(Job.Nifti,2);
end
if nargin<4
	LogName = NaN;
end
dd_initcnode(Job)

% Wrap around to get the Masks for all Series
for n = 1:numel(SeriesNr)
	[Mask{n} Target{n}] = getmeanmask(Job, SubjNr, SeriesNr(n), LogName);
end


function [Mask Target] = getmeanmask(Job, SubjNr, SeriesNr, LogName)
%%
% This function gets the mask for a single series

if iscell(Job.BETOpts.Str)
	BETOpts	= Job.BETOpts.Str{SubjNr,SeriesNr};
	BETMeth = Job.BETMenu.Str{Job.BETMenu.Val(SubjNr,SeriesNr)};
else
	BETOpts	= Job.BETOpts.Str;
	BETMeth = Job.BETMenu.Str{Job.BETMenu.Val};
end
Realign = ~strcmp(Job.RealignMenu.Str{Job.RealignMenu.Val}, 'none');

%-- Just get the existing mask and do nothing else
[b0s DWIs] = getb0imgs(Job, SubjNr, SeriesNr);
if isempty(LogName) || isnan(LogName(1))
	[Mask Target] = dd_getmeanmask(b0s, Realign, BETOpts, LogName, SeriesNr);
	return
end

%-- Create or get the b0/DWI-mask
switch BETMeth
	case 'mean(b0)'
		
		%-- Get mean-b0 volumes
		[Mask Target] = dd_getmeanmask(b0s, Realign, BETOpts, LogName, SeriesNr);
		
	case 'mean(DWI)'
		
		% Strategy: replace the b0-based mask with DWI-based one. In this way
		% the mask always has the same name, which is useful for interactive bet-usage

		%-- First get existing DWI or new mean-b0 volumes (without printing)
		[Mask Target Brain New] = dd_getmeanmask(b0s, Realign, BETOpts, NaN, SeriesNr);

		if New	% We only have b0-based mask -> replace it with DWI-based one

			%-- Make a working copy of the first 20 DWIs (i.e. enough to get a good enough SNR for BET)
			TmpDir = tempname;
			mkdir(TmpDir)
			for n = 1:min(size(DWIs,1),20)
				[~, Name Ext] = fileparts(DWIs(n,:));
				DWIs_temp{n}  = fullfile(TmpDir, [Name Ext]);
				copyfile(spm_file(DWIs(n,:),'ext','.*'), TmpDir)
			end
			DWVols = spm_vol(char(DWIs_temp));

			%-- Realign DWI -> mean-b0 and get a coregistered mean-DWI volume
			spm_realign(DWVols, struct('sep',2, 'graphics',0));		% Store the realignment info in the temporary header files
			spm_reslice(DWVols, struct('which',0, 'mean',1))		% Save a mean-DWI-volume for coregistration
			x = spm_coreg(spm_file(DWVols(1).fname,'prefix','mean'), Target, struct('graphics',0));	% x = from mean DWI to mean b0-images
			for n = 1:numel(DWVols)									% Update the realignment info in the header
				RM = spm_matrix(x) * spm_get_space(DWVols(n).fname);
				spm_get_space(DWVols(n).fname, RM);
			end
			spm_reslice([spm_vol(Target); DWVols], struct('which',1, 'mean',0, 'prefix',''))	% Reslice the DWI-volumes in the space of the mean b0

			%-- Make and print a new DWI mask and copy it into Mask & recreate Brain
			[DWIMask MeanDWI] = dd_getmeanmask(char(DWVols.fname), false, BETOpts, LogName, SeriesNr);
			MaskVol			  = spm_read_vols(spm_vol(DWIMask));
			spm_write_vol(spm_vol(Mask), MaskVol);
			spm_write_vol(spm_vol(Brain), spm_read_vols(spm_vol(Target)) .* MaskVol);

			%-- Save the mean DWI image and clean up
			movefile(spm_file(MeanDWI,'ext','.*'), fileparts(DWIs(1,:)), 'f')
			rmdir(TmpDir, 's')

		else	% Print the pre-existing DWI mask
			
			dd_getmeanmask(b0s, Realign, BETOpts, LogName, SeriesNr);
			
		end
		
	otherwise
		
		error('Unknown Option: %s', BETMeth)
		
end
if isempty(Mask)
	return
end

% Save info to the logfile
MaskVol = spm_read_vols(spm_vol(Mask));
VoxVol	= abs(det(spm_get_space(Mask)))/1000;				% Voxel volume in ml
FIDLog	= fopen([LogName(1:end-2) 'txt'], 'a');
fprintf(FIDLog, 'S%d\t%s-mask:\tVolume (ml) =\t%g\n', SeriesNr, BETMeth, sum(MaskVol(:)) * VoxVol);
fclose(FIDLog);


function [b0Imgs DWImgs] = getb0imgs(Job, SubjNr, SeriesNr)
%%
% This function gets the b0-images for a single series
%
% Input:  dd_basicproc Job
% Output: filelist (SPM-style)

if isempty(Job.Nifti(SubjNr,SeriesNr).Files)
    warning('DIDI:OpenFile', ['Cannot find any images in: ' Job.Nifti(SubjNr,SeriesNr).Path])
	[b0Imgs DWImgs] = deal([]);
	return
end
for n = 1:numel(Job.Nifti(SubjNr,SeriesNr).Files)
	FList{n} = fullfile(Job.Nifti(SubjNr,SeriesNr).Path, Job.Nifti(SubjNr,SeriesNr).Files{n});
    b(n)	 = getfield(dti_get_dtidata(FList{n}), 'b');
end
b0Sel  = b<=50;   % ==0;
b0Imgs = char(FList{b0Sel});
DWImgs = char(FList{~b0Sel});

if isempty(b0Imgs)
    warning('DIDI:OpenFile', ['Cannot find any b0-images in: ' Job.Nifti(SubjNr,SeriesNr).Path])
end
