function [x output] = peunwarp_fminunc_regularised(targetVol, sourceVol, order, maskVol)

% FUNCTION [x output] = peunwarp_fminunc_regularised(targetFile, sourceFile, order, maskFile)

global bX bY dbY bZ target8 source mask krn xscale regstrength XFlip

%-- Parse the input
if nargin < 3
    order = [9 9 9];
end
if nargin<3
	maskVol = '';
end

%-- Constants
fwhm		= 2;				% Histogram smoothing
regstrength = single(.1);		%
xscale      = single(1e3);

%-- Read the images from disk
if isempty(maskVol)
	mask = ones(sourceVol.dim, 'uint8');
else
	mask = uint8(spm_read_vols(spm_vol(maskVol)));
end
targetVol = spm_vol(targetVol);
sourceVol = spm_vol(sourceVol);
target	  = spm_read_vols(targetVol);
source	  = spm_read_vols(sourceVol);
if targetVol.mat(1)<0
	disp('LR flipped image:')
	disp(targetVol.mat)
	XFlip = true;
else
	XFlip = false;
end

%-- Ceil out susceptibility artefacts in the source image intensities to improve binning
[Dum SortInd] = sort(source(:));
CeilFrac = round(0.999 * numel(source));	% Ceil the source intensity at the 99th percentile
% mask(SortInd(CeilFrac:end)) = 0;			% This is not optimal, it would be better to mask the source instead of the target
source(SortInd(CeilFrac:end)) = source(SortInd(CeilFrac));

%-- Ceil the target image intensities to improve binning
[Dum SortInd] = sort(target(:));
CeilFrac = round(0.999 * numel(target));	% Ceil the target intensity at the 99th percentile
% mask(SortInd(CeilFrac:end)) = 0;
target(SortInd(CeilFrac:end)) = target(SortInd(CeilFrac));

%-- Scale images to 0-255 and cast to the required class-type
targetMin = min(target(:));
targetMax = max(target(:));
target8	  = uint8(round((target - targetMin) * 255 / (targetMax - targetMin)));
sourceMin = min(source(:));
sourceMax = max(source(:));
source	  = single((source - sourceMin) * 255 / (sourceMax - sourceMin));

%-- Create the dct basis functions
bX  = single(spm_dctmtx(sourceVol.dim(1), order(1)));
bY  = single(spm_dctmtx(sourceVol.dim(2), order(2)));
dbY = single(spm_dctmtx(sourceVol.dim(2), order(2), 'diff'));
bZ  = single(spm_dctmtx(sourceVol.dim(3), order(3)));
bX  = bX / bX(1, 1);
dbY = dbY / bY(1, 1);
bY  = bY / bY(1, 1);
bZ  = bZ / bZ(1, 1);

%-- From spm_coreg
% lim = ceil(2 * fwhm);
lim = ceil(fwhm);
krn = smoothing_kernel(fwhm, -lim:lim);
krn = single(krn / sum(krn));

%-- Find the optimal unwarping parameters
x0 = zeros(prod(order), 1);					% Start with zero deformation
% x = fminsearch(@peunwarp_fminunc_costwrapper, x0, optimset('Display', 'iter-detailed',  'TolFun', 1e-10));
% x = pewarp_spm_powell(x, xi, tolsc, maxiter, 'mex_pewarpcost', bX, bY, dbY, bZ, target8, source, mask, krn, xscale);
% if verLessThan('optim','4.2')
% 	DisplayOpt = 'iter';
% else
% 	DisplayOpt = 'iter-detailed';
% end
% x = fminunc(@peunwarp_fminunc_regularised_costwrapper, x0, optimset('Display', DisplayOpt, ...
%     'OutputFcn', @peunwarp_fminunc_regularised_plot, 'LargeScale', 'off', 'Diagnostics', 'on', ...
%     'FinDiffType', 'central', 'MaxFunEvals', 1000000, 'TolX', 1e-10));
if exist('minFunc','file')~=2
	addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'minFunc_2012')))
	RemPath = true;
else
	RemPath = false;
end
% [x, ~, ~, output] = minFunc(@peunwarp_fminunc_regularised_costwrapper, x0, ...
% 	struct('Display','iter', 'outputFcn',@peunwarp_fminunc_regularised_plot, ...
% 		   'numDiff',1, 'MaxIter',200, 'MaxFunEvals',1000000, 'TolX',1e-12));
[x dum dum output] = minFunc(@peunwarp_fminunc_regularised_costwrapper, x0, ...
	struct('outputFcn',@peunwarp_fminunc_regularised_plot, ...
		   'numDiff',1, 'MaxIter',50, 'MaxFunEvals',1000000, 'progTol',1e-12));
if RemPath
	rmpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'minFunc_2012')))
end

%% END %%


function f = peunwarp_fminunc_regularised_costwrapper(x)

global bX bY dbY bZ target8 source mask krn xscale regstrength

f = mex_pewarpcost_regularised(single(x), bX, bY, dbY, bZ, target8, source, mask, krn, xscale, regstrength);


function Stop = peunwarp_fminunc_regularised_plot(x, Type, Iteration, varargin)

% Old required arguments (v2009):	outputFcn(x,infoStruct,state,varargin{:})
% New required arguments (v2012):	outputFcn(x,iterationType,i,funEvals,f,t,gtd,g,d,optCond,varargin{:});
%
% Changed:	.._plot(x, Info, State, varargin)	=> .._plot(x, Type, Iteration, varargin)
%			Info.iteration						=> Iteration
%			State								=> Type

Stop = false;

if Iteration==0 && strcmp(Type, 'iter')
	return;				% fminunc calls twice (no need to plot again)
end

global bX bY dbY bZ target8 source mask krn xscale regstrength XFlip

[Cost JHist Unwarped Def] = mex_pewarpcost_regularised(single(x), bX, bY, dbY, bZ, target8, source, mask, krn, xscale, regstrength);

h = spm_figure('FindWin', 'Interactive');
if ishandle(h)
	
	if Iteration==0
		spm_figure('Clear', h);
		spm('FigName', 'PE-unwarp optimization', h);
	end
	set(0, 'CurrentFigure', h)

	PrevCost = get(findobj(h, 'Tag', 'MIOpt'), 'YData');
	plot(0:length(PrevCost), [PrevCost -Cost], 'Tag', 'MIOpt')
	set(gca, 'XTick', unique(round(get(gca,'XTick'))))
	xlabel('Iteration')
	ylabel('Mutual information')
	
end	

h = spm_figure('FindWin', 'Graphics');
if ishandle(h)
	
	MeanMaxDef = [mean(abs(Def(logical(mask)))) max(abs(Def(logical(mask))))];
	setappdata(h, 'MeanMaxDef', MeanMaxDef);

	XGrid = {'XTickLabel',[],'YTickLabel',[], 'XGrid','On', 'GridLineStyle',':', 'XColor','c'};
	XRange = size(target8,1):-1:1;
	if XFlip
		XRange = fliplr(XRange);
	end
	set(0, 'CurrentFigure', h)
	set(h, 'Visible', 'On')
	spm_figure('Clear', h);
	spm('FigName', 'PE-unwarp results', h);
	colormap('gray')
	
	subplot(2, 2, 1)
	imagesc([shiftdim(target8(floor(end/2),:,:))'; target8(XRange,:,floor(end/2))])
	set(gca, XGrid{:})
	title('Target')
	axis xy image
	
	subplot(2, 2, 2)
	MaxSrc = max(source(floor(end/2),:));
	imagesc([shiftdim(source(floor(end/2),:,:))'; source(XRange,:,floor(end/2))], [0 MaxSrc])
	set(gca, XGrid{:})
	title('Source')
	axis xy image
	
	subplot(2, 2, 3)
	imagesc([shiftdim(Unwarped(floor(end/2),:,:))'; Unwarped(XRange,:,floor(end/2))], [0 MaxSrc])
	set(gca, XGrid{:})
	title(['Unwarped (' Type ' ' num2str(Iteration) ')'])
	axis xy image
	
	subplot(2, 2, 4)
	DiffIm = single(Unwarped) - source;
	imagesc([shiftdim(DiffIm(floor(end/2),:,:))'; DiffIm(XRange,:,floor(end/2))], [-MaxSrc MaxSrc])
	set(gca, XGrid{:})
	title('Unwarped - Source')
	axis xy image
	
%	DefIm = shiftdim(Def(floor(end/2), :, :))';
% 	subplot(3, 2, 5)
% 	imagesc(-DefIm, [-1 1] * max(abs([DefIm(:); eps])))
% 	set(gca, XGrid{:})
% 	title('Deformation in PE-direction')
% 	axis xy image
% 	text(0.04, 0.8, sprintf('mean %.1f\nmax %.1f', MeanMaxDef), 'Units','Normalized')
% 	try colorbar('SouthOutside'), end
% 	
% 	subplot(3, 2, 6)
% 	% imagesc(JHist, [0 50])
% 	imagesc(log(1+JHist))
% 	title('Joint histogram')
% 	set(gca, 'XTick',[], 'YTick',[])
% 	xlabel('Target')
% 	ylabel('Source')
% 	axis xy image
	
	refresh(h)

end


% From spm_coreg
function krn = smoothing_kernel(fwhm,x)

% Variance from FWHM
s = (fwhm/sqrt(8*log(2)))^2+eps;

% The simple way to do it. Not good for small FWHM
% krn = (1/sqrt(2*pi*s))*exp(-(x.^2)/(2*s));

% For smoothing images, one should really convolve a Gaussian
% with a sinc function.  For smoothing histograms, the
% kernel should be a Gaussian convolved with the histogram
% basis function used. This function returns a Gaussian
% convolved with a triangular (1st degree B-spline) basis
% function.

% Gaussian convolved with 0th degree B-spline
% int(exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= -0.5..0.5)
% w1  = 1/sqrt(2*s);
% krn = 0.5*(erf(w1*(x+0.5))-erf(w1*(x-0.5)));

% Gaussian convolved with 1st degree B-spline
%  int((1-t)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= 0..1)
% +int((t+1)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t=-1..0)
w1  =  0.5*sqrt(2/s);
w2  = -0.5/s;
w3  = sqrt(s/2/pi);
krn = 0.5*(erf(w1*(x+1)).*(x+1) + erf(w1*(x-1)).*(x-1) - 2*erf(w1*x   ).* x)...
      +w3*(exp(w2*(x+1).^2)     + exp(w2*(x-1).^2)     - 2*exp(w2*x.^2));

krn(krn<0) = 0;