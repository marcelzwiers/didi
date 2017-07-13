function [PWPar FEVal] = peunwarp_estimate(Ref, Src, Order, Mask)

% function [PWPar FEVal] = peunwarp_estimate(Ref, Src, Order, [Mask])
%
% Estimates the deformation field (warping) parameters that map the source
% (e.g. DWVol) images onto the (resliced space of the) reference image (e.g. T1Vol)
%
% Ref, Src and Mask must already be coregistered.

if nargin<3 || isempty(Order)
    % Order = [10 16 10];
	Order = [9 9 9];
end
if nargin<4
	Mask = '';
end
if isstruct(Ref), Ref = Ref.fname; end
if isstruct(Src), Src = Src.fname; end

% Load the reference image into memory and create a temporary copy header
RefHdr		  = spm_vol(Ref);
RefVol		  = spm_read_vols(RefHdr);
SrcHdr		  = spm_vol(Src);
Ext			  = ['.' spm_file(Ref,'ext')];
rRef		  = [tempname Ext];			% Ref-copy in local tempdir is fast & easy
rRefHdr		  = rmfield(SrcHdr,'pinfo');
rRefHdr.fname = rRef;

% Filter the temporary copy before reslicing to avoid aliasing artefacts (i.e. decimate)
R	   = RefHdr.mat \ SrcHdr.mat;		% Maps voxel coordinates of Src onto Ref
SampR  = max(abs(R(1:3,1:3)), 2);		% Voxelsize Src in reference space
Win	   = bsxfun(@times, fwin(RefHdr.dim(1),SampR(1)) * fwin(RefHdr.dim(2),SampR(2))', ...
						shiftdim(fwin(RefHdr.dim(3),SampR(3)), -2));
RefVol = ifftn(Win .* fftn(RefVol), 'symmetric');	% Avoid potential tiny imaginary roundoff errors

% Reslice the filtered copy onto the source image using trilinear interpolation (same as in pewarp)
%spm_reslice(char(Src, rRef), struct('mean', 0, 'which', 1))
% NB: the code below is probably correct but spm_reslice('interp',1) also use trilinear interpolation (see spm_cfg_realign)
rRefVol = zeros(SrcHdr.dim);
for n = 1:SrcHdr.dim(3)
	rRefVol(:,:,n) = spm_slice_vol(RefVol, R * spm_matrix([0 0 n]), SrcHdr.dim(1:2), 1);
end
rRefHdr = spm_write_vol(rRefHdr, rRefVol);

% Reslice the mask onto the source image using nearest-neighbour interpolation (must remain logical)
if ~isempty(Mask)
	rMask			= [tempname Ext];			% Ref-copy in local tempdir is fast & easy
	rMaskHdr		= rmfield(SrcHdr,'pinfo');
	rMaskHdr.fname	= rMask;
	rMaskHdr.dt(1)	= spm_type('uint8');
	MaskVol			= spm_read_vols(spm_vol(Mask));
	rMaskVol		= zeros(MaskHdr.dim);
	for n = 1:MaskHdr.dim(3)
		rMaskVol(:,:,n) = spm_slice_vol(MaskVol, R * spm_matrix([0 0 n]), SrcHdr.dim(1:2), 0);
	end
	Mask = spm_write_vol(rMaskHdr, rMaskVol);
end

% Estimate the unwarping field
tic
[PWPar Output] = peunwarp_fminunc_regularised(rRefHdr, SrcHdr, Order, Mask);
FEVal.Elapsed  = toc;
FEVal.NrEvals  = Output.funcCount;

% Save the unwarping parameters
[SrcPath SrcName] = fileparts(Src);
save(fullfile(SrcPath, [SrcName '_pwparam.mat']), 'PWPar')


%% END

function w = fwin(nn, SR)

% w = 1D zero-padded Tukey filter window matched to the sampling ratio
% 
% Filter window is defined in three sections: taper, constant, taper
% Period of the taper is defined as 1/2 period of a sine wave.

if SR < 1.25
	w = ones(nn,1);			% No filtering if we do no (or hardly any) downsampling
	return
end

% Constants
per = 0.25;					% = r/2, r=0.5 (= 1/4 taper + 1/2 constant + 1/4 taper)

m = round(nn/SR);			% Size of resampling (Src) window (must be < nn)
if rem(m,2)==0;
	n = m + 1;				% Make sure the Tukey window is symmetric, i.e. has an odd length 'n'
else
	n = m;
end

% Create the Tukey filter window of length n
t  = linspace(0,1,n)';
tl = floor(per*(n-1))+1;
th = n-tl+1;
tw = [((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];

% Zero-pad the filter window to the full (length nn) k-space
w  = ifftshift([zeros(ceil((nn-n)/2), 1); tw; zeros(floor((nn-n)/2), 1)]);
