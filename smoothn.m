function sVol = smoothn(Vol, FWHM)

% FUNCTION sVol = smoothn(Vol, FWHM)
%
% SMOOTHN uses fft to perform Gaussian smoothing
%
% INPUT
%	Vol  - n-dimensional data array that needs to be smoothed
%	FWHM - Full-width-half-max of the Gaussian convolution kernel expressed
%	       in voxel indices. FWHM must be a single scalar or an array with
%	       the same length as the number of dimensions of Vol
%
% OUTPUT
%	sVol - Smoothed n-dimensional data array
%
% Marcel, 26-11-2010

if ndims(Vol)>1 && numel(FWHM)==1
	FWHM(1:ndims(Vol)) = FWHM;
end

Win = fwin(size(Vol,1), FWHM(1));
for n = 2:ndims(Vol)
	Win = bsxfun(@times, Win, shiftdim(fwin(size(Vol,n), FWHM(n)), -(n-1)));
end
if isreal(Vol)
	sVol = ifftn(Win .* fftn(Vol), 'symmetric');	% Avoid potential tiny imaginary roundoff errors
else
	sVol = ifftn(Win .* fftn(Vol));
end


function Win = fwin(n, FWHM)

if FWHM==0
	Win = ones(n,1);
	return
end

x	= [-floor(n/2):floor(n/2)]';
a	= sqrt(log(2));				% 0.5 = exp(-a^2), x = FWHM/2
Win = exp(-(2*a*x/FWHM).^2);
Win = fft(ifftshift(Win));		% Transform the window to the frequency domain
Win = Win / Win(1);				% Make sure we don't loose energy
if rem(n,2)==0
	Win(round(end/2)) = [];		% Make the number of matrix elements even again
end
