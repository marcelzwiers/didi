function [RBTM APar D2TM] = dd_basicproc_realign(Job, SubjNr, RBTM, LogName, Constr, Tuning, HWS)

% DD_BASICPROC_REALIGN is a key internal function of dd_basicproc that is made
% available externally to allow distributed computing. Create a DWGrandTVol,
% coregister it with the T1 and rotate all DWTVols to this coregistered grand
% target. NB: This function does change volume headers but not the diffusion
% gradient directions (this happens in DD_BASICPROC_WARP).
%
% USAGE
%	[RBTM APar D2TM] = DD_BASICPROC_REALIGN(Job, SubjNr, RBTM, LogName, Constr, Tuning, HWS)
%
% It returns the rigid-body / affine transformation matrix relative to b0-space
% RBTM (i.e. the realignment parameters, not the motion parameters!) and/or,
% depending on the Job, the affine transformation parameters (for eddy-current
% correction) and the transformation to undistorted space.
%
% Marcel, 21-10-2014
%
% See also: DD_BASICPROC, DD_BASICPROC_REALIGNPAR

% Defaults
if nargin < 7
	HWS = [];
end

% Initialize main outcome parameters
NrSeries	= numel(Job.Nifti(SubjNr,:));
[APar D2TM] = deal(repmat({},1,NrSeries));

% Return if the pipeline failed at an earlier stage
if any(cellfun(@isempty,{Job.Nifti(SubjNr,:).Files}))
	return
end

% Launch a graphics window if we do a distributed job
dd_initcnode(Job)

%--- Realign DWTVols to create a DWGrandTVol. NB: realignment info is stored in memory in DWTVol
[~, Meanb0] = dd_basicproc_getmask(Job, SubjNr);		% Get the b0-target(s) for each Series
DWTVol		= spm_vol(char(Meanb0));					% SPM targets for each Series
DWTVol		= spm_realign(DWTVol, struct('quality',1, 'sep',2, 'rtm',1));	% This should not change DWTVol on disk...

%--- Compute the transformation (R) from the DWGrandTVol to the T1 reference image
T1Img = char(Job.Nifti(SubjNr,1).T1Text);
if Job.T1Box.Val && ~isempty(T1Img)
	
	if NrSeries>1
		spm_reslice(DWTVol, struct('which',0))			% Write-out a temporary grand mean target volume
		[Pth Nme Ext] = fileparts(DWTVol(1).fname);
		DWGrandTVol	  = spm_vol(fullfile(Pth, ['mean' Nme Ext]));
	else
		DWGrandTVol	  = DWTVol;
	end
	x = spm_coreg(spm_vol(T1Img), DWGrandTVol, struct('graphics',true));
	if NrSeries>1
		delete(fullfile(Pth, ['mean' Nme '.*']))		% Clean-up the temporary grand mean target
	end
	R = inv(spm_matrix(x));								% Transformation from DWGrandTVol to T1Img
	x = spm_imatrix(R);
	myspm_print(LogName, 'S*: Coregistration DWI -> T1', Job.Nifti(SubjNr).Path)
	FIDLog = fopen([LogName(1:end-2) 'txt'], 'a');
	fprintf(FIDLog, 'S*\tCoregistration DWI->T1:\tRBTPar =\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n', ...
					x(1:3), x(4:6)*180/pi);
	fclose(FIDLog);
	
	TDVal = logical([Job.MDBox.Val Job.NormDBox.Val Job.NormABox.Val Job.FABox.Val Job.RABox.Val Job.ModeBox.Val]);
	if Job.NormBox.Val && any(TDVal) && NrSeries==1		% Realign T1 with T1-template to make spm_preproc more robust
		NTempl	= fullfile(spm('Dir'), 'canonical', 'avg152T1.nii');
		MT1		= spm_get_space(T1Img);					% The original T1-orientation
		PT1		= spm_realign(char(NTempl, T1Img), struct('graphics',0,'rtm',0));
		R		= (PT1(2).mat / MT1) * R;				% Add T1=>T1-template rotation to R
		spm_get_space(T1Img, PT1(2).mat);				% Rotate T1=>T1-template
	end
	
else
	
	R = eye(4);											% Transformation from DWGrandTVol to T1Img
	
end

%--- Compute the rigid-body / affine transformation parameters for each Series
for SeriesNr = 1:NrSeries
	[RBTM{SeriesNr} APar{SeriesNr} D2TM{SeriesNr}] = realignseries(DWTVol(SeriesNr), R, ...
		Job, SubjNr, SeriesNr, RBTM{SeriesNr}, LogName, Constr, Tuning, HWS);
end


function [RBTM APar D2TM] = realignseries(DWTVol, R, Job, SubjNr, SeriesNr, RBTM, LogName, Constr, Tuning, HWS)

[APar D2TM] = deal([]);

% Get the DWVol headers, the DW info and the b0/DWI-mask
DWVol = spm_vol(makeflist(Job.Nifti(SubjNr,SeriesNr)));	% SPM source input
for n = 1:numel(DWVol)
	D(n) = dti_get_dtidata(DWVol(n).fname);
end

%--> Compute the transformation (R) from DWTVol -> DWGrandTVol -> (MNI realigned) T1 image
M	   = spm_get_space(DWTVol.fname);	% DWTVol is in register on disk with the raw DWI/meanb0 target volume
R	   = R * DWTVol.mat / M;			% DWTVol is in register in memory with DWGrandTVol
DWTVol = spm_vol(DWTVol.fname);			% Restore DWTVol now that the info is stored in R

%--> Compute and set the transformation (RBTM) from DWVol -> DWTVol in all volume headers
if ~strcmp(Job.RealignMenu.Str{Job.RealignMenu.Val}, 'none')
	
	% Realign the b0-images to DWTVol (Mask volumes should already be in register with DWTVol)
	b0s = [D.b]<50;
	spm_realign([DWTVol; DWVol(b0s)], struct('quality',1, 'sep',2, 'rtm',0));	% spm_realign should do just fine here
	
	% Compute the transformation from DWVol -> DWTVol / mean-DWI
	[RBTM APar] = dd_basicproc_realignpar(Job, SubjNr, SeriesNr, LogName, Constr, Tuning, HWS, RBTM);
	
	% Compute the transformation from mean-DWI -> DWTVol and set all headers -> DWTVol
	s	   = find(~b0s);										% = the index of DWVol(~b0Sel); numel(s)==size(RBTM,3)
	FIDLog = fopen([LogName(1:end-2) 'txt'], 'a');
	switch Job.RealignMenu.Str{Job.RealignMenu.Val}
		
		case 'sq DT-error'										% D2TM = DWVol -> DWGrandTVol = WIP & not tested
			% Use mean undistorted image to find transform (D2TM) between b0/ref- and diffusion-space
			fprintf('\nSampling the mean distorted space')
			df_reslice(APar, struct('mask',0, 'hold',1, 'mean',1, 'which',2),[],@df_unwarp_graphwip);	% Create a mean undistorted diffusion-weighted image
			fprintf('\n--> Coregistering the mean distorted space\n')
			DWMean		= spm_file(APar.P(1).fname, 'prefix','mean');
			D2TPar		= spm_coreg(DWMean, DWTVol, struct('params',spm_imatrix(eye(4)), 'graphics',true));		% NB: DWTVol in memory -> DWGrandTVol
			D2TM		= spm_matrix(D2TPar);
			D2TPar(4:6) = D2TPar(4:6)*180/pi;
			delete([DWMean(1:end-3) '*'])						% Clean up the tmp mean-DWI-volume
			myspm_print(LogName, ['S' num2str(SeriesNr) ': Coregistration DWI -> b0'], Job.Nifti(SubjNr).Path)
			fprintf(FIDLog, 'S%g\tCoregistration DWI->b0:\tRBTPar =\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n', ...
							SeriesNr, D2TPar(1:6));
			disp(['DW-space to b0 (mm & deg):' sprintf('\t%.2f', D2TPar)])
			
			[S LWarn] = mywarning('Off');						% Get rid of QFORM0 warning from affine transformations
			for n = 1:numel(s)									% Move outside switch block?
				RM = RBTM(:,:,n) * DWVol(s(n)).mat;				% Not sure, check!
				spm_get_space(DWVol(s(n)).fname, RM);			% Store the realignment info in the header: This may not be 100% correct because APar needs to be applied first?
			end
			mywarning(S, LWarn)
			
		case 'mutual info'										% RBTM = DWVol -> DWTVol
			[S LWarn] = mywarning('Off');						% Get rid of QFORM0 warning from affine transformations
			for n = 1:numel(s)
				RM = RBTM(:,:,n) * DWVol(s(n)).mat;
				spm_get_space(DWVol(s(n)).fname, RM);			% Store the realignment info in the header on disk
			end
			mywarning(S, LWarn)
			
		case 'sq error'											% RBTM = DWVol -> mean-DWI
			% Use mean-DWI to find transform between b0/ref- and diffusion-space
			MDWVol = spm_vol(spm_file(DWVol(s(1)).fname, 'prefix','mean'));
			x	   = spm_coreg(MDWVol, DWTVol.fname, struct('graphics',true));	% x = from mean DWI to mean b0-images
			myspm_print(LogName, ['S' num2str(SeriesNr) ': Coregistration DWI -> b0'], Job.Nifti(SubjNr).Path)
			fprintf(FIDLog, 'S%g\tCoregistration DWI->b0:\tRBTPar =\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n', ...
							SeriesNr, x(1:3), x(4:6)*180/pi);
			disp(['DW-space to b0 (mm & deg):' sprintf('\t%.2f', x(1:3), x(4:6)*180/pi)])
			
			[S LWarn] = mywarning('Off');						% Get rid of QFORM0 warning from affine transformations
			for n = 1:numel(s)
				RM = spm_matrix(x) * RBTM(:,:,n) * DWVol(s(n)).mat;
				spm_get_space(DWVol(s(n)).fname, RM);			% Store the realignment info in the header
			end
			mywarning(S, LWarn)
			
	end
	DWVol = spm_vol(char(DWVol.fname));							% Update DWVol.mat
	
	% Plot and log the results (does not work for sq DT-error)
	RBTPar_all = zeros(numel(DWVol), 12);
	for n = 1:numel(DWVol)
		RBTPar_all(n,:) = spm_imatrix(DWVol(n).mat / M);		% NB: realignment parameters = inv(motion parameters)
	end
	dd_basicproc_plotparameters(RBTPar_all, DWVol, LogName, ['S' num2str(SeriesNr) ': Realignment DWI (patched)'], Job.ParallelBox.Val)
	myset(HWS, 'Visible', 'Off')
	fprintf(FIDLog, 'S%d\tRealignment (patched):\trange(RBTPar) =\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tRMS_FD(RBTPar) =\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
					SeriesNr, myrange(RBTPar_all(:,1:6)), myrms(RBTPar_all(:,1:6)));
	fclose(FIDLog);
	save(spm_file(DWTVol.fname, 'prefix','rp_', 'ext','.txt'), 'RBTPar_all', '-ascii');
	
end

%--> Left-append R and update all headers (DWTVol,DWVol,Mask,Brain) -> DWGrandTVol
if Job.T1Box.Val
	[S LWarn] = mywarning('Off');					% Get rid of QFORM0 warning from affine transformations
	spm_get_space(DWTVol.fname, R*M);
	for n = 1:numel(DWVol)
		RM = R * DWVol(n).mat;						% Rotation from DWVol voxel to DWGrandTVol world
		spm_get_space(DWVol(n).fname, RM);			% Map DWVol to DWGrandTVol world
		if strfind(DWVol(n).fname, 'artc_')			% Idem for the original (raw) files
			spm_get_space(strrep(DWVol(n).fname,'artc_',''), RM);
		end
	end
	Mask = dd_basicproc_getmask(Job, SubjNr, SeriesNr);	% Don't print the mask again
	spm_get_space(char(Mask), R*M);						% _brain_mask
	spm_get_space(strrep(char(Mask),'_mask',''), R*M);	% _brain
	mywarning(S, LWarn)
end


%% ------------------- Auxillary Functions -------------------


function Range = myrange(X)

Range      = max(X) - min(X);
Range(4:6) = Range(4:6)*180/pi;


function Displace = mydiff(X)

Displace      = sum(abs(diff(X)));
Displace(4:6) = Displace(4:6)*180/pi;


function RMS = myrms(X)

% Gives more weight to spike movements compared to mydiff
X(:,4:6) = X(:,4:6)*180/pi;
RMS      = sqrt(mean(diff(X).^2));


function FList = makeflist(Data)
%
% Input:  Job.Data (.Path & .Files)
% Output: filelist (SPM-style)

for n = 1:numel(Data.Files)
    FList{n} = fullfile(Data.Path, Data.Files{n});
end
FList = char(FList);


function myset(varargin)

% Robust against closed figures
if ishandle(varargin{1})
    set(varargin{:})
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
