function [pwVol Data] = peunwarp_write(Vol, PWPar, Order, TMat, Data, APar, D2TM)

% FUNCTION [Vol Data] = peunwarp_write(Vol, PWPar, Order, TMat, Data)
%
% Unwarps the series of source images in the PE-direction (dim(2))
%
% INPUT:
%	Vol		- A character array or structure containing header information (see spm_vol)
%	PWPar	- Unwarping parameters (see peunwarp_estimate) that map the DWTVol (b0) image
%			  onto the (realigned and resliced) T1 reference image
%	Order	- Order of spatial basis functions (same as used in peunwarp_estimate)
%	TMat	- Additional affine transformation matrices of size(4,4,numel(Vol)) that
%			  map the DWI images onto the DWTVol image
%	Data	- A xyzt-matrix that is unwarped (i.e. passed to spm_bsplinc) instead of
%			  the Vol data. Useful for mapping associated (e.g. weighting) volumes
%
% OUTPUT:
%	Vol		- A vector of structures containing header information of the
%			  unwarped images (prefix 'pw')
%	Data	- Unwarped xyzt data
%
% Marcel, 6-8-2010

% PRIVATE INPUT (see DD_BASICPROC_REALIGNPAR & DD_BASICPROC):
%	APar	- Slicewise eddy-current & head motion correction parameters (from df_unwarp)
%	D2TM	- Rigid body transformation matrix from mean diffusion space to b0 (target) space

% Defaults
Deg	   = 4 * [1 1 1];	% Degree of B-spline (from 0 to 7) along different dimensions
PreFix = 'pw';
if nargin==0
	pwVol = PreFix;
	return
end

% Check the input
Vol = spm_vol(Vol);

% Construct the deformation field
Bx  = spm_dctmtx(Vol(1).dim(1), Order(1));
By  = spm_dctmtx(Vol(1).dim(2), Order(2));
Bz  = spm_dctmtx(Vol(1).dim(3), Order(3));
Bx  = Bx / Bx(1, 1);
By  = By / By(1, 1);
Bz	= Bz / Bz(1, 1);
Def = reshape(spm_get_def(Bx, By, Bz, PWPar*1e3), Vol(1).dim);

% Construct the sampling grid
[X Y Z] = ndgrid(1:Vol(1).dim(1), 1:Vol(1).dim(2), 1:Vol(1).dim(3));

spm_progress_bar('Init', numel(Vol), 'PE-unwarping', 'Images to write')
pwVol = Vol;
for n = 1:numel(Vol)

	%--> Add extra transformation(s). NB: inverse transform of coordinates = transform of image
	XYZ = [X(:) Y(:) Z(:) ones(numel(X),1)]';
	if nargin>3
		if (nargin==4 || nargin==5 || (nargin>5 && isempty(APar))) && ~isempty(TMat)	% Use TMat
			XYZ = (Vol(n).mat \ TMat(:,:,n) * Vol(n).mat) \ XYZ;						% Transform TMat to voxelspace
		elseif nargin==7 && ~isempty(APar)					% Use APar (from df_unwarp) and D2TM
			Sel = [1:6+3*Vol(n).dim(3)] + (n-1)*(6+3*Vol(n).dim(3));
			XYZ = df_transf_coord(Vol(n), APar.beta(Sel), XYZ, [], D2TM, false);
		end
	end
	
	%--> Resample and save the data to disk
	if nargin>4 && ~isempty(Data)
		% Data(:,:,:,n) = reshape(spm_sample_vol(Data(:,:,:,n), XYZ(1,:), XYZ(2,:)+Def(:)', XYZ(3,:), 1), Vol(n).dim);
		pwData		  = spm_bsplins(spm_bsplinc(Data(:,:,:,n),Deg), XYZ(1,:), XYZ(2,:)+Def(:)', XYZ(3,:), Deg);
		pwData(~isfinite(pwData)) = 0;		% spm_bsplins gives NaN's for out-of-volume sampling points
		Data(:,:,:,n) = reshape(pwData, Vol(n).dim);	% Overwrite 'Data' with the result
		pwVol		  = [];					% We don't save the result to disk
	else
		% pwData	   = spm_sample_vol(Vol(n), XYZ(1,:), XYZ(2,:)+Def(:)', XYZ(3,:), 1);
		pwData		   = spm_bsplins(spm_bsplinc(Vol(n),Deg), XYZ(1,:), XYZ(2,:)+Def(:)', XYZ(3,:), Deg);
		pwData(~isfinite(pwData)) = 0;		% spm_bsplins gives NaN's for out-of-volume sampling points
		pwVol(n).fname = spm_file(Vol(n).fname, 'prefix', PreFix);
	    %pwVol(n).mat   = Vol(1).mat;		% Update the orientation info in the header
		pwVol(n).pinfo = [0;0;0];
		[S LWarn]	   = mywarning('Off');	% Get rid of QFOMR0 warnings
		pwVol(n)	   = spm_write_vol(pwVol(n), reshape(pwData, Vol(n).dim));
		mywarning(S, LWarn)
	end
	
	% Update the waitbar
	spm_progress_bar('Set',n)
	
end