function DNFiles = dd_basicproc_denoising(Files, Method, Rician, LogName, SeriesNr)

% DD_BASICPROC_DENOISING is a key internal function of dd_basicproc that
% is made available externally to allow distributed computing.
%
% FUNCTION DNFiles = DD_BASICPROC_DENOISING(FNames, Method, Rician, LogName)
%
% Marcel, 17-3-2014
%
% See also: DD_BASICPROC, MAINDWIDENOISING

% Defaults
if nargin<2 || isempty(Method)
	Method = 'LPCA';
end
if nargin<3 || isempty(Rician)
	Rician = false;
end
if nargin<4 || isempty(LogName)
	Verbose = false;			% Don't display the figure
else
	Verbose = true;				% Display the figure
end
if nargin<5 || isempty(SeriesNr)
	SeriesNr = 1;
end
if ~exist('DWIDenoisingLPCA','file')
	addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'DWIDenoisingPackage_r01_pcode')))
end
Beta	  = 1;					% Default smoothing parameter (not for AONLM and LPCA)
PRadius   = 1;					% Default patch-radius
[S LWarn] = mywarning('Off','MATLAB:maxNumCompThreads:Deprecated');
NThreads  = maxNumCompThreads;	% Ignore Intel's hyperthreading capabilities
mywarning(S, LWarn)

% Make sure Files is a SPM-style charlist and not a cell-in-cell
if iscell(Files) && numel(Files)==1
	Files = Files{1};
end
Files = char(Files);
if isempty(Files)
	DNFiles = '';
	return
end

% Read the data
Hdrs = spm_vol(Files);
% Vol = spm_read_vols(Hdrs);					% Not robust against realigned images
Vol	 = zeros([Hdrs(1).dim numel(Hdrs)]);
for n = 1:numel(Hdrs)
	Vol(:,:,:,n) = spm_read_vols(Hdrs(n));
	D(n) = dti_get_dtidata(Hdrs(n).fname);
end

% Scale, denoise and unscale the data
Range = [min(Vol(:)) max(Vol(:))];
Vol	  = 255 * (Vol - Range(1)) / diff(Range);	% Scale the data to [0,255]
if ~any(strcmp(Method, {'AONLM' 'LPCA'}))
	Dirs   = vertcat(D.g);						% Matrix containing DTI gradient directions in world-space
	HFinal = DWINoiseEstimation(Vol, Dirs, Rician, Verbose);
	if isnan(HFinal)
		warning('DIDI:Denoise:NoiseEstimation', 'Error during noise estimation')
	end
end
switch Method
	case 'ONLM'
		DNVol = DWIDenoisingORNLM(Vol, HFinal, Beta, PRadius, 5, Rician, NThreads, Verbose);
	case 'AONLM'
		DNVol = DWIDenoisingAONLM(Vol, PRadius, 3, Beta , Rician, NThreads, Verbose);
	case 'Multires-ONLM'
		DNVol = DWIDenoisingORNLMMultires(Vol, Dirs, HFinal, Beta, PRadius, 3, Rician, NThreads, Verbose);
	case 'ODCT'
		DNVol = DWIDenoisingODCT(Vol, HFinal, Beta, Rician, NThreads, Verbose);
	case 'PRINLM'
		DNVol = DWIDenoisingPRINLM(Vol, HFinal, Beta, Rician, NThreads, Verbose);
	case 'LPCA'
		DNVol = DWIDenoisingLPCA(Vol, Beta, Rician, NThreads, Verbose);
	otherwise
		error('Unknown denoising method: %s', Method)
end
close
Vol	  = diff(Range) *   Vol / 255 + Range(1);		% Scale the data back to the original range
DNVol = diff(Range) * DNVol / 255 + Range(1);		% Scale the data back to the original range
DNVol(DNVol<0 | isnan(DNVol)) = 0;					% I'm not sure if this is necessary / applicable

% Save the denoised data
DNHdrs = rmfield(Hdrs, 'pinfo');
for n = 1:numel(DNHdrs)
	DNHdrs(n).fname	= spm_file(DNHdrs(n).fname, 'prefix',[Method '_']);
	spm_write_vol(DNHdrs(n), DNVol(:,:,:,n));
	dti_get_dtidata(DNHdrs(n).fname, D(n));			% Update/create the dw-info
end
DNFiles	= char(DNHdrs.fname);

% Print the results
if Verbose
	
	[S LWarn] = mywarning('Off', 'SPM:noDisplay');
	HG		  = spm_figure('GetWin', 'Graphics');
	mywarning(S, LWarn)
	spm_figure('Clear', HG)
	
	b0Sel = [D.b]<=50;								% ==0;
	b0s	  = find(b0Sel);
	DWIs  = find(~b0Sel);
	zi	  = round(size(Vol,3)/2);
	PSize = [0.27 0.27];
	Pos	  = [0.08 0.38 0.68];
	if Hdrs(1).mat(1)<0
		Vol	  = flipdim(Vol,1);
		DNVol = flipdim(DNVol,1);
	end

	H		 = axes('Position', [Pos(1) Pos(3) PSize], 'Visible','Off', 'Parent',HG);
	RawSlice = Vol(:, :, zi, b0s(1));
	Range	 = [min(RawSlice(:)) max(RawSlice(:))];
	imagesc(flipud(RawSlice'), 'Parent',H, Range)
	axis(H, 'image')
	set(H, 'XTick',[], 'YTick',[])
	ylabel(H, sprintf('g=[%.2f  %.2f  %.2f]',D(b0s(1)).g))
	title(H, 'Raw')

	H		 = axes('Position', [Pos(2) Pos(3) PSize], 'Visible','Off', 'Parent',HG);
	DNSlice  = DNVol(:, :, zi, b0s(1));
	imagesc(flipud(DNSlice'), 'Parent',H, Range)
	axis(H, 'image', 'off')
	title(H, sprintf('Denoised (%s)',Method))

	H		 = axes('Position', [Pos(3) Pos(3) PSize], 'Visible','Off', 'Parent',HG);
	imagesc(flipud((RawSlice-DNSlice)'), 'Parent',H)
	axis(H, 'image', 'off')
	title(H, 'Noise')

	H		 = axes('Position', [Pos(1) Pos(2) PSize], 'Visible','Off', 'Parent',HG);
	RawSlice = Vol(:, :, zi, DWIs(1));
	Range	 = [min(RawSlice(:)) max(RawSlice(:))];
	imagesc(flipud(RawSlice'), 'Parent',H, Range)
	axis(H, 'image')
	set(H, 'XTick',[], 'YTick',[])
	ylabel(H, sprintf('g=[%.2f  %.2f  %.2f]',D(DWIs(1)).g))

	H		 = axes('Position', [Pos(2) Pos(2) PSize], 'Visible','Off', 'Parent',HG);
	DNSlice  = DNVol(:, :, zi, DWIs(1));
	imagesc(flipud(DNSlice'), 'Parent',H, Range)
	axis(H, 'image', 'off')

	H		 = axes('Position', [Pos(3) Pos(2) PSize], 'Visible','Off', 'Parent',HG);
	imagesc(flipud((RawSlice-DNSlice)'), 'Parent',H)
	axis(H, 'image', 'off')

	H		 = axes('Position', [Pos(1) Pos(1) PSize], 'Visible','Off', 'Parent',HG);
	RawSlice = mean(Vol(:, :, zi, :), 4);
	Range	 = [min(RawSlice(:)) max(RawSlice(:))];
	imagesc(flipud(RawSlice'), 'Parent',H, Range)
	axis(H, 'image')
	set(H, 'XTick',[], 'YTick',[])
	ylabel(H, 'mean(DWI)')

	H		 = axes('Position', [Pos(2) Pos(1) PSize], 'Visible','Off', 'Parent',HG);
	DNSlice  = mean(DNVol(:, :, zi, :), 4);
	imagesc(flipud(DNSlice'), 'Parent',H, Range)
	axis(H, 'image', 'off')

	H		 = axes('Position', [Pos(3) Pos(1) PSize], 'Visible','Off', 'Parent',HG);
	NSlice	 = std(Vol(:, :, zi, :) - DNVol(:, :, zi, :), [], 4);
	imagesc(flipud(NSlice'), 'Parent',H)
	axis(H, 'image')
	set(H, 'XTick',[], 'YTick',[], 'YAxisLocation','Right')
	ylabel(H, 'std(Noise)', 'Rotation',-90, 'VerticalAlignment','Bottom')

	H		 = axes('Position', [0.1 0.95 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0, 0.5, sprintf('S%d: Denoising (%s) ../%s', SeriesNr, Method, spm_file(Hdrs(1).fname,'basename')), 'Parent',H)
	H		 = axes('Position', [0.05 0.00 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0, 0.5, datestr(now), 'Parent',H)
	colormap(HG, 'gray')
	
	[S LWarn] = mywarning('Off','spm:spm_jobman:NotInitialised');
	spm_print(LogName)
	mywarning(S, LWarn)
	
end
