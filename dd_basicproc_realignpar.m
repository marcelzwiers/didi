function [RBTM APar] = dd_basicproc_realignpar(Job, SubjNr, SeriesNr, LogName, Constr, Tuning, HWS, RBTM)

% DD_BASICPROC_REALIGNPAR is an internal function of dd_basicproc that is made
% available externally to allow distributed computing. This function does not
% change anything and only returns the rigid-body transformation matrices RBTM
% (i.e. the realignment parameters, not the motion parameters!) and the
% skull-stripped brain mask.
%
% Marcel, 26-8-2010
%
% See also: DD_BASICPROC, DD_BASICPROC_REALIGN, DD_BASICPROC_REALIGNWARP

% Return if the pipeline failed at an earlier stage
APar = [];
if isempty(Job.Nifti(SubjNr,SeriesNr).Files) || strcmp(Job.RealignMenu.Str{Job.RealignMenu.Val}, 'none')
	RBTM = [];
	return
end

if nargin < 7
	HWS = [];
end
if nargin > 7
	Stage = 'patched';
else
	Stage = 'raw';
end

% Launch a graphics window if we do a distributed job
dd_initcnode(Job)

% Get the mask and DTI info for all Series
for n = 1:numel(Job.Nifti(SubjNr,SeriesNr).Files)
	D(n) = orderfields(dti_get_dtidata(fullfile(Job.Nifti(SubjNr,SeriesNr).Path, Job.Nifti(SubjNr,SeriesNr).Files{n})));
end
b0Sel		  = [D.b]<=50;   % ==0;
[Mask Meanb0] = dd_basicproc_getmask(Job, SubjNr, SeriesNr);	% Don't print the mask again

% Compute RBTM if realignment has to be performed
fprintf('\n-> %s\n', Job.Nifti(SubjNr,SeriesNr).Path)

DWVol  = spm_vol(makeflist(Job.Nifti(SubjNr,SeriesNr)));
DWVol  = DWVol(~b0Sel);						% Realign the DW images only
DWTVol = spm_vol(char(Meanb0));				% = SPM target input
FIDLog = fopen([LogName(1:end-2) 'tsv'], 'a');
RBTPar = repmat(spm_imatrix(eye(4)), size(DWVol));
if ~exist('RBTM','var') || isempty(RBTM)
	RBTM	  = repmat(eye(4), [1 1 numel(DWVol)]);
	EmptyRBTM = true;
else
	EmptyRBTM = false;
end
switch Job.RealignMenu.Str{Job.RealignMenu.Val}
	
	case 'sq DT-error'
		disp('-> Performing realignment and eddy current correction...')
		if strcmp(Stage,'raw')													% PATCH does not support slicewise transformations
			Constr.sptl.order = 1;												% Use 1st order spatial basis functions
		end
		Tuning.imask.vol  = spm_read_vols(spm_vol(char(Mask)));					% Add a brain-mask
		[D DimFlip]		  = dd_rotategradients(D(~b0Sel), inv(DWTVol.mat));		% Recompute the DW info in voxel coordinates (=> df_unwarp does not consider the volume orientation)
		BVecs			  = diag(DimFlip) * vertcat(D.g)';						% This flip is because Volkmar didn't do his signs right (see dd_rotategradients)
		APar = df_unwarp(DWVol, BVecs', Constr, Tuning, @df_unwarp_graphwip);	% Give the gradient vectors in voxel space
		myspm_print(LogName, sprintf('S%d: Realignment (%s DT-error)',SeriesNr, Stage), ...
					Job.Nifti(SubjNr,SeriesNr).Path)
		fprintf(FIDLog, 'S%g\tRealignment (%s DT-error):\tNIter =\t%d\tMSS =\t%g\n', SeriesNr, Stage, ...
						length(APar.mss), APar.mss(end));
		NSlc = DWVol(1).dim(3);
		for n = 1:numel(DWVol)
			if strcmp(Stage,'raw')
				Betas	= APar.beta([1:6+3*NSlc] + (n-1)*(6+3*NSlc))';	% Not supported by dd_patch so take the mean over the volume
				Betas	= [Betas(1:6) 1 1 1 0 0 0] + ...
						  [0 mean(Betas(6 + [1:NSlc] + 2*NSlc)) 0  0 0 0 ...
						   0 mean(Betas(6 + [1:NSlc] + NSlc)) 0 ...
						   0 mean(Betas(6 + [1:NSlc])) 0];		% = Affine transformation in voxel coordinates
			else
				Betas	= APar.beta([1:6] + (n-1)*(6+3*NSlc))';
			end
			RBTPar(n,:) = spm_imatrix(inv(spm_matrix(Betas)));	% Reverse the inverse (sampling) transformation of APar.beta (see df_transf_coord)
		end

	case 'mutual info'
		myset(HWS, 'Visible', 'On');
		for n = 1:numel(DWVol)
			fprintf('\n-> coregistering volume: %g\n', n);
			if strcmp(Stage,'raw') || EmptyRBTM || ~isempty(strfind(DWVol(n).fname, [filesep 'artc_']))
				RBTPar(n,:) = spm_coreg(DWVol(n), DWTVol, struct('params',spm_imatrix(RBTM(:,:,n)), 'graphics',false));
			else												% We already have computed the uncorrected (raw) coregistration, no need to do it again
				RBTPar(n,:) = spm_imatrix(RBTM(:,:,n));
			end
			mywaitbar(n/numel(DWVol), HWS, 'Coregistering DW-volumes...')
			fprintf('Realignment parameters: %s\n', sprintf('\t%f', RBTPar(n,:)))
		end

	case 'sq error'												% NB: two-pass affine realignment wrt mean-DWI (thus not DWTVol)
		myset(HWS, 'Visible', 'On');
		for n = 2:numel(DWVol)
			RBTM(:,:,n)	 = spm_affreg(DWVol(1), DWVol(n), struct('regtype','rigid'), RBTM(:,:,n));
			DWVol(n).mat = RBTM(:,:,n) * DWVol(n).mat;			% Temporarily apply the realignment for creation of a mean-DWI volume
			mywaitbar(n/numel(DWVol), HWS, 'Realigning DW-volumes (pass 1)...')
		end
		spm_reslice(DWVol, struct('which',0, 'mean',1))			% Save a mean-DWI-volume
		MDWVol = spm_vol(spm_file(DWVol(1).fname, 'prefix','mean'));
		disp('Realignment:')
		for n = 1:numel(DWVol)
			DWVol(n).mat = RBTM(:,:,n) \ DWVol(n).mat;			% Undo the temporary realignment
			RBTPar(n,:)	 = spm_imatrix(spm_affreg(MDWVol, DWVol(n), struct('regtype','rigid'), RBTM(:,:,n)));
			mywaitbar(n/numel(DWVol), HWS, 'Realigning DW-volumes (pass 2)...')
			fprintf('[%.3f%s]\n', RBTPar(n,1), sprintf('\t%.3f', RBTPar(n,2:end)))
		end

	otherwise
		error('Realignment option not recognized')
		
end
myset(HWS, 'Visible', 'Off')

% Create the relative transformation matrices (x y z pitch rol yaw scale[x y z] shear[x y z])
for n = 1:size(RBTPar,1)
	RBTM(:,:,n) = spm_matrix(RBTPar(n,:));
end

% Print the result
if strcmp(Stage,'raw')
	dd_basicproc_plotparameters(RBTPar, DWVol, LogName, sprintf('S%g: Realignment DWI (%s)',SeriesNr,Stage), Job.ParallelBox.Val)
	fprintf(FIDLog, 'S%g\tRealignment (%s DWI):\trange(RBTPar) =\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tRMS_FD(RBTPar) =\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
					SeriesNr, Stage, myrange(RBTPar(:,1:6)), myrms(RBTPar(:,1:6)));
end
fclose(FIDLog);


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


function FList = makeflist(DSet)
%
% Input:  Job.DSet (.Path & .Files)
% Output: filelist (SPM-style)

for n = 1:numel(DSet.Files)
    FList{n} = fullfile(DSet.Path, DSet.Files{n});
end
FList = char(FList);


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
