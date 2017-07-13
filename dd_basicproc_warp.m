function [TgtImgs WImgs] = dd_basicproc_warp(Job, SubjNr, LogName, APar, D2TM)

% DD_BASICPROC_WARP is a key internal function of dd_basicproc that
% is made available externally to allow distributed computing.
%
% USAGE
%	[TgtImgs WImgs] = DD_BASICPROC_WARP(Job, SubjNr, LogName, APar, D2TM)
%
% Marcel, 10-11-2014
%
% See also: DD_BASICPROC, DD_BASICPROC_REALIGN

dd_initcnode(Job)
PWPar = [];

for SeriesNr = 1:size(Job.Nifti,2)
	[TgtImgs{SeriesNr} WImgs{SeriesNr} PWPar] = warp(Job, SubjNr, SeriesNr, LogName, ...
													 APar{SeriesNr}, D2TM{SeriesNr}, PWPar);
end


function [TgtImgs WImgs PWPar] = warp(Job, SubjNr, SeriesNr, LogName, APar, D2TM, PWPar)

% Return if the pipeline failed at an earlier stage
if isempty(Job.Nifti(SubjNr,SeriesNr).Files)
	TgtImgs = [];
	WImgs	= [];
	return
end

% Launch a graphics window if we do a distributed job
HG = spm_figure('FindWin', 'Graphics');

% Get some DWVol info and a b0-mask
T1Img		  = Job.Nifti(SubjNr,SeriesNr).T1Text;
[Mask Meanb0] = dd_basicproc_getmask(Job, SubjNr, SeriesNr);
Mask		  = char(Mask);
DWTVol		  = spm_vol(char(Meanb0));								% = SPM target input
DWVol		  = spm_vol(makeflist(Job.Nifti(SubjNr,SeriesNr)));		% = SPM source input
RBTM		  =	zeros(4,4,numel(DWVol));
for n = 1:numel(DWVol)
	RBTM(:,:,n) = DWVol(n).mat / DWTVol.mat;
	D(n)		= dti_get_dtidata(DWVol(n).fname);
end
b0Sel = [D.b]<=50;			% =b0;

%--> Compute the total weight (= WLog * WVox * WSlc)
if strcmp(Job.EstMenu.Str{Job.EstMenu.Val}, 'PATCH')
	disp('Loading & constructing total weights...')
	load(fullfile(Job.Nifti(SubjNr,SeriesNr).Path, 'PATCH.mat'))
	WTot = bsxfun(@times, PATCH.WVox, PATCH.WSlc);
	if PATCH.WLog
		WTot = PATCH.WLog .* WTot;
	end
	clear PATCH
else
	WTot = [];
end

%--> Estimate and write realigned PE-unwarped images or reslice the realigned images, i.e. warp(DWVol,WTot,Brain,Mask,DWTVol)
if Job.T1Box.Val && Job.PEUnwarpBox.Val && ~isempty(T1Img)
	
	%-- Estimate the co-registered T1 and mean b0-image
	disp('Estimate PE-unwarping...')
	Order = str2num(Job.PEUnwarpText.Str);
	if ~Job.SinglePEUnwarpBox.Val || SeriesNr==1
		[PWPar FEval] = peunwarp_estimate(T1Img, DWTVol, Order);	% Unwarp the meanb0-image to the T1 reference image. Use mask (e.g. to elimate the influence of shifted fat-rims)?
		myspm_print(LogName, ['S' num2str(SeriesNr) ': Unwarping b0 -> T1'], ...
					sprintf('Nr evaluations: %d, elapsed time: %d', FEval.NrEvals, FEval.Elapsed))
	end
	% Get the absolute in-mask deformation and save this as a QA-parameter
	if ishandle(HG)
		FIDLog = fopen([LogName(1:end-2) 'txt'], 'a');
		fprintf(FIDLog, 'S%g\tabs(PEUnwarpDef):\tmean =\t%g\tmax =\t%g\n', SeriesNr, getappdata(HG,'MeanMaxDef'));
		fclose(FIDLog);
	end
	
	%-- Unwarp DWVol + WTot
	PreFix = peunwarp_write;							% = 'pw'
	peunwarp_write(DWVol( b0Sel), PWPar, Order, RBTM(:,:, b0Sel));
	peunwarp_write(DWVol(~b0Sel), PWPar, Order, RBTM(:,:,~b0Sel), [], APar, D2TM);
	if ~isempty(WTot)
		WTot	   = permute(WTot, [4 3 2 1]);			% The weights need to be in xyzt-order
		[Dum WTot] = peunwarp_write(DWVol(~b0Sel), PWPar, Order, RBTM(:,:,~b0Sel), WTot, APar, D2TM);
		WTot	   = ipermute(WTot, [4 3 2 1]);			% Restore the tzyx-order
	end
	
	%-- Unwarp Brain + Mask + DWTVol to prevent realignment and brain extraction when calling [..]_getmeanmask in a later stage
	peunwarp_write(strrep(Mask,'_mask',''), PWPar, Order); % Write unwarped Brain, but name = pwmean[..] instead of meanpw[..]
	peunwarp_write(Mask, PWPar, Order);					% Write unwarped Mask, but name = pwmean[..] instead of meanpw[..] and not binary yet
	writemeanmask(DWVol(b0Sel), Mask, PreFix)			% Write DWTVol (= MeanImg) and rename and binarize Brain and Mask
	
elseif ~strcmp(Job.RealignMenu.Str{Job.RealignMenu.Val}, 'none')
	
	if ~isempty(WTot)									% Create a sampling grid for WTot
		Dim		= DWTVol.dim;
		[Z Y X] = ndgrid(1:Dim(3), 1:Dim(2), 1:Dim(1));	% Z is the fastest varying dimension
		XYZ		= [X(:) Y(:) Z(:) ones(prod(Dim(1:3)),1)]';
	end
	if ~isempty(APar)									% RealignMenu = 'sq DT-error' = success

		%-- Unwarp DWVol + WTot
		PreFix = 'u';
		disp('Saving unwarped & resliced DW volumes to disk...')
		spm_reslice(DWVol(b0Sel), struct('which',2, 'mean',0, 'prefix',PreFix))
		df_reslice(APar, struct('mask',0, 'hold',1, 'mean',0, 'which',1), [], @df_unwarp_graphwip, D2TM);
		if ~isempty(WTot)
			disp('Saving resliced weight volumes to disk...')
			spm_progress_bar('Init', size(WTot,1), 'Reslicing weights', 'Volumes to reslice')
			for n = 1:size(WTot,1)
				Betas	  = APar.beta([1:6+3*Dim(3)] + (n-1)*(6+3*Dim(3)));
				tXYZ	  = df_transf_coord(DWTVol, Betas, XYZ, [], D2TM); % Undistorted space + extra transformation
				WTot(n,:) = spm_sample_vol(shiftdim(WTot(n,:,:,:)), tXYZ(3,:), tXYZ(2,:), tXYZ(1,:), 1);
				spm_progress_bar('Set',n)
			end
		end
		
	else
		
		%-- Reslice DW images + WTot
		PreFix = 'r';									% Use the default value
		disp('Saving resliced DW volumes to disk...')
		spm_reslice(DWVol, struct('which',2, 'mean',0))
		if ~isempty(WTot)
			spm_progress_bar('Init', size(WTot,1), 'Reslicing weights', 'Volumes to reslice')
			for n = 1:size(WTot,1)
				tXYZ	  = (DWTVol.mat \ RBTM(:,:,n) * DWTVol.mat) \ XYZ;
				WTot(n,:) = spm_sample_vol(shiftdim(WTot(n,:,:,:)), tXYZ(3,:), tXYZ(2,:), tXYZ(1,:), 1);
				spm_progress_bar('Set',n)
			end
		end
	end
	clear X Y Z XYZ tXYZ
	
	%-- Reslice Brain + Mask + DWTVol to prevent realignment and brain extraction when calling [..]_getmeanmask in a later stage
	spm_reslice(strrep(Mask,'_mask',''), struct('which',2, 'mean',0, 'prefix',PreFix))	% Reslice Brain, but name = premean[..] instead of meanpre[..]
	spm_reslice(Mask, struct('which',2, 'mean',0, 'prefix',PreFix))						% Reslice Mask, but name = premean[..] instead of meanpre[..] and not binary yet
	writemeanmask(DWVol(b0Sel), Mask, PreFix)											% Write DWTVol (= MeanImg) and rename and binarize Brain and Mask
	
elseif strcmp(Job.RealignMenu.Str{Job.RealignMenu.Val}, 'none')
	
	PreFix = '';
	
end

% Ready with spatial alignments, create SPM Diffusion toolbox variables (Volkmar Glauche)
TgtImgs = prepend({DWVol(:).fname}', PreFix);			% = DTI-tbx target input

%--> Copy the DW info to the resliced images and (re)set spatial orientations
rD = D;													% A rotated copy of D (belonging to TgtImgs)
if PreFix
	for n = 1:numel(TgtImgs)
		rD(n)	  = dd_rotategradients(D(n), RBTM(:,:,n));	% Rotate the gradient vector
		rD(n).mat = DWTVol.mat;							% This is our new reference space
		spm_get_space(fullfile(Job.Nifti(SubjNr,SeriesNr).Path,Job.Nifti(SubjNr,SeriesNr).Files{n}), D(n).mat);	% Reset the orientation of the original input images
		spm_get_space(TgtImgs{n}, DWTVol.mat);			% Make sure all targets have the same orientation
		dti_get_dtidata(TgtImgs{n}, rD(n));				% Update/create the dw-info
	end
end

%--> Save total weight to nifti-files for use in dd_wregress
if ~isempty(WTot)
	s = find(~b0Sel);
	for n = 1:size(WTot,1)
		WVolHdr			= rmfield(spm_vol(TgtImgs{s(n)}), 'pinfo');
		WVolHdr.fname	= postpend(WVolHdr.fname, '_wt');
		WVolHdr.dt(1)	= spm_type('uint8');
		WVolHdr.descrip	= 'Voxelwise total weight (for robust diffusion model estimation)';
		WVolHdr			= spm_write_vol(WVolHdr, permute(WTot(n,:,:,:), [4 3 2 1]));	% The weights need to be in xyzt-order
		WImgs{n}		= WVolHdr.fname;
	end
else
	WImgs = '';
end
WImgs = char(WImgs);


%% ------------------- Auxillary Functions -------------------


function writemeanmask(DWVol, Mask, PreFix)
%
% Write mean* and give proper name and datatype to *_brain and *_brain_mask to
% prevent realignment and brain extraction when calling [..]_getmeanmask

% Parse the input
PreFixVol  = prepend({DWVol.fname}',PreFix);
PreFixMask = prepend(Mask, PreFix);
Mask	   = strrep(PreFixMask, [PreFix 'mean'], ['mean' PreFix]);	% Our desired unwarped Mask

% Write MeanImg
spm_reslice(PreFixVol, struct('mean',1, 'which',0))

% Rename *_brain and *_brain_mask
movefile(strrep(PreFixMask,'_mask',''), strrep(Mask,'_mask',''), 'f')	% premean[..]_brain		 -> meanpre[..]_brain
movefile(PreFixMask, Mask, 'f')											% premean[..]_brain_mask -> meanpre[..]_brain_mask
if strcmp(spm_file(Mask,'ext'),'img')
	movefile(spm_file(strrep(PreFixMask,'_mask',''),'ext','hdr'), spm_file(strrep(Mask,'_mask',''),'ext','hdr'), 'f')
	movefile(spm_file(PreFixMask,'ext','hdr'), spm_file(Mask,'ext','hdr'), 'f')
end

% Change the datatype of the Mask
MaskHdr		  = spm_vol(Mask);
MaskHdr.dt(1) = spm_type('int16');						% Same as in BET
spm_write_vol(MaskHdr, round(spm_read_vols(MaskHdr)));	% Make the unwarped/resliced mask binary again


function FList = makeflist(DSet)
%
% Input:  Job.DSet (.Path & .Files)
% Output: filelist (SPM-style)

for n = 1:numel(DSet.Files)
    FList{n} = fullfile(DSet.Path, DSet.Files{n});
end
FList = char(FList);


function PO = prepend(PI, Pre)
%
% Input: filelist (SPM-style) or cellarray

for n = 1:size(PI,1)
    [pth,nm,xt] = fileparts(char(PI(n,:)));
    PO{n}       = fullfile(pth, [Pre nm xt]);
end
if ischar(PI)
	PO = char(PO);
end


function PO = postpend(PI, Post)
%
% Input: filelist (SPM-style) or cellarray

for n = 1:size(PI,1)
    [pth,nm,xt] = fileparts(char(PI(n,:)));
    PO{n,1}     = fullfile(pth, [nm Post xt]);
end
if ischar(PI)
	PO = char(PO);
end


function myspm_print(LogName, HdrTxt, FtrTxt)

% Robust against closed figures
HG = spm_figure('FindWin', 'Graphics');
if ishandle(HG)
	HD = axes('Position', [0.05 0.95 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0.75, 0.5, datestr(now), 'Parent',HD)
	text(0, 0.5, HdrTxt, 'Parent',HD)
	HD = axes('Position', [0.05 0 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0, 0.5, FtrTxt, 'Parent',HD)
	spm_print(LogName)
end
