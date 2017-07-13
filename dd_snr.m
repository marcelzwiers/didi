function [SNR BGStd GhStd Spikes tSNR] = dd_snr(Vol, Mask, FNames, Title)

% DD_SNR is a helper function that computes SNR levels for DD_BASICPROC
%
% FUNCTION [SNR BGStd GhStd Spikes tSNR] = DD_SNR(Vol, Mask, FNames, Title)
%
% INPUT
%	Vol		- Filenames of the xyz-data volumes or tzyx data array
%	Mask	- Filename of the xyz-mask volume or zyx logical array. If left
%			  empty a simple intensity mask is constructed
%	FNames	- Filenames of the xyz-data volumes
%	Title	- Title of one of the axes (if not-empty a graphical plot is made)
%
% OUTPUT
%	SNR		- Estimated SNR from background noise: [SNR_SIG SNR_SD]
%	BGStd	- SD of background noise per volume per slice
%	GhStd	- SD of ghost-image per volume
%	Spikes	- Nr of slices with abnormal high background signal (i.e. from RF spikes)
%	tSNR	- mean temporal SNR in mask
%
% Marcel, 27-01-2011
%
% See also: DD_BASICPROC

if nargin<4
	Title = '';
end
if nargin<3
	FNames = '';
end
if nargin<2
	Mask = [];
end
if isempty(FNames) && ischar(Vol)
	FNames = Vol;
end

% Load data if necessary
if ischar(Mask) || iscellstr(Mask)
	Mask = logical(permute(spm_read_vols(spm_vol(char(Mask))), [3 2 1]));
end
if ischar(Vol)
	% Vol = permute(spm_read_vols(spm_vol(Vol)), [4 3 2 1]); % Not robust against realigned images
	Hdr  = spm_vol(Vol(1,:));
	Vols = zeros([Hdr.dim size(Vol,1)]);
	for n = 1:size(Vol,1)
		Vols(:,:,:,n) = spm_read_vols(spm_vol(Vol(n,:)));
	end
	Vol = permute(Vols, [4 3 2 1]);						% tzyx-volume
	clear Vols
end
if isempty(Mask)										% Use a simple threshold mask
	MeanVol = smoothn(shiftdim(mean(Vol,1)), 5);
	Mask	= MeanVol > 0.3*max(MeanVol(:));
end
R		= 2;											% Create a spherical kernel with radius R
Kernel	= ones(2*R+1, 2*R+1, 2*R+1);
[I J K] = ind2sub(size(Kernel), find(Kernel));
Kernel(sqrt((I-R-1).^2 + (J-R-1).^2 + (K-R-1).^2) > R) = 0;
Air		= bsxfun(@times, Vol, shiftdim(spm_erode(1-Mask,Kernel),-1));	% Zeros will be omitted in BGStd & GhStd

% Construct a background mask
d	= 10;
Box	= [1 d size(Vol,3)-d+1 size(Vol,3)
	   1 d size(Vol,4)-d+1 size(Vol,4)];				% The outer corner boxes
BG0	= Air(:, :, [Box(1,1):Box(1,2) Box(1,3):Box(1,4)], ...
				[Box(2,1):Box(2,2) Box(2,3):Box(2,4)]);	% NB: Slices can be zero after reslicement

% Construct a ghost-image mask
[Z Y X] = ind2sub(size(Mask), find(Mask));
Slb	= [1 d size(Vol,3)-d+1 size(Vol,3)
	   min(Z) max(Z) min(X) max(X)];					% A square slab in front+behind the brain
BG1	= Air(:, Slb(2,1):Slb(2,2), [Slb(1,1):Slb(1,2) Slb(1,3):Slb(1,4)], ...
			 Slb(2,3):Slb(2,4));

% Compute the background mean and std per volume per slice and the ghost-STD per volume
BGAvg = zeros(size(BG0,1), size(BG0,2));
BGStd = zeros(size(BG0,1), size(BG0,2));
GhStd = zeros(size(BG0,1));
for ti = 1:size(BG0,1)
	for zi = 1:size(BG0,2)
		DSel = BG0(ti,zi,:)>=1;							% Omit the arteficial zeros from Siemens and reslicement
		if sum(DSel) > 1
			BGAvg(ti,zi) = mean(BG0(ti,zi,DSel));
			BGStd(ti,zi) = std(BG0(ti,zi,DSel));
		end
	end
	GhStd(ti) = std(BG1(ti, BG1(ti,:)>0));				% NB: Slices can be zero after reslicement
end
if any(BGStd(:)==0)
	fprintf('Background noise voxels were all zero in %g/%g slices (likely due to realignment, than this warning is unimportant)\n', ...
			sum(BGStd(:)==0), numel(BGStd))
	BGAvg(BGAvg==0) = median(BGAvg(:));					% Replace them with a rough estimate
	BGStd(BGStd==0) = median(BGStd(:));					% Replace them with a rough estimate
end

% Estimate the SNR: Henkelman RM, 1985. Measurement of signal intensities in the presence of noise in MRI images. Med Phys 12:232-33
%SNR = mean(mean(Vol(:,Mask))) * [1.253/mean(BG0(BG0>=1)) 0.655/mean(BGStd(:))];
SNR = mean(mean(Vol(:,Mask))) * [1.253/mean(BGAvg(:)) 0.655/mean(BGStd(:))];
if nargout>4
	tSNR = mean(mean(Vol(:,Mask)) ./ std(Vol(:,Mask)));
end

% Robustly detect spikes in the background
Thresh = median(BGAvg(:)) + 4 * 1.4826 * median(abs(BGAvg(:)-median(BGAvg(:))));	% Thresh = mean + 4 SD (MAD estimated)
Spikes = [sum(BGAvg(:) > Thresh) numel(BGAvg)];

% Display the results
HG = spm_figure('FindWin', 'Graphics');
if ~isempty(HG) && ~isempty(Title)
	
	spm_figure('Clear', HG)
	set(HG, 'visible', 'on')
	[Nr Int] = hist(reshape(BG0, numel(BGStd), [])', 25);	% The histograms of the raw background data
	
    ax = axes('Position',[0.1 0.77 0.8 0.16],'Parent',HG,'Visible','off');
	y = 1;
	text(0,y, 'Estimation of background noise','FontSize',16,'FontWeight','Bold','Parent',ax)
	y = y - 0.13;
	text(0,y, [fileparts(FNames(1,:)) ':'],'FontSize',9,'Interpreter','none','Parent',ax)
	y = y - 0.09;
	for i = 1:min([size(FNames,1) 10])
		[Dum Name Ext] = fileparts(FNames(i,:));
		if length(Name)>45
			Name = [Name(1:5) '[...]' Name(end-35:end)];
		end
		text(0,y, sprintf('%-4.0f%s%s',i,Name,Ext),'FontSize',9,'Interpreter','none','Parent',ax)
		y = y - 0.09;
	end
	if size(FNames,1) > 10
		text(0,y,'................ etc','FontSize',9,'Parent',ax)
	end
	
	ax = axes('Position',[0.63 0.75 0.3 0.17],'Parent',HG,'Visible','off');
	StdVol = shiftdim(max(std(Vol,[],1),[],2));
	Mask   = double(shiftdim(max(Mask)));
	if exist('Hdr','var') && Hdr.mat(1)<0
		StdVol = fliplr(StdVol);
		Mask   = fliplr(Mask);
	end
	StdVol = 1.5 * StdVol / max(StdVol(:));					% Make the image a bit brighter [0:1.5]+0.3
	StdVol = repmat(StdVol, [1 1 3]) + ...
			 bsxfun(@times, Mask, shiftdim([0.3;0;0],-2));	% mask=red, alpha=0.3
	StdVol(StdVol>1) = 1;
	image(StdVol, 'Parent', ax)
	axis(ax, 'xy', 'image')
	hold(ax, 'on')
	Box = Box + [-0.5 0.5 -0.5 0.5
				 -0.5 0.5 -0.5 0.5];
	plot(ax, [Box(2,1) Box(2,1) Box(2,3) Box(2,3)
			  Box(2,1) Box(2,1) Box(2,3) Box(2,3)
			  Box(2,2) Box(2,2) Box(2,4) Box(2,4)
			  Box(2,2) Box(2,2) Box(2,4) Box(2,4)
			  Box(2,1) Box(2,1) Box(2,3) Box(2,3)], ...
			 [Box(1,1) Box(1,3) Box(1,1) Box(1,3)
			  Box(1,2) Box(1,4) Box(1,2) Box(1,4)
			  Box(1,2) Box(1,4) Box(1,2) Box(1,4)
			  Box(1,1) Box(1,3) Box(1,1) Box(1,3)
			  Box(1,1) Box(1,3) Box(1,1) Box(1,3)], 'c-')
	Slb = Slb + [-0.5 0.5 -0.5 0.5
				 -0.5 0.5 -0.5 0.5];
	plot(ax, [Slb(2,3) Slb(2,3)
			  Slb(2,3) Slb(2,3)
			  Slb(2,4) Slb(2,4)
			  Slb(2,4) Slb(2,4)
			  Slb(2,3) Slb(2,3)], ...
			 [Slb(1,1) Slb(1,3)
			  Slb(1,2) Slb(1,4)
			  Slb(1,2) Slb(1,4)
			  Slb(1,1) Slb(1,3)
			  Slb(1,1) Slb(1,3)], 'r-')
	
	ax = subplot(4,1,2, 'Parent',HG);
	plot(ax, Int, Nr)
	XLim2 = double(max(Int(:)));
	xlim(ax, [0 XLim2])
	YLim2 = 1.1 * max(max(Nr(2:end,:)));
	ylim(ax, [0 YLim2])										% The first bin contains the zeros
	text(0.75*XLim2, 0.8*YLim2, {sprintf('BG   = %.1f (%.1f)', mean(BG0(BG0>1)), mean(BGStd(:)))
								 sprintf('SNR = %.1f (%.1f)', SNR)
								 sprintf('GNR = %.1f', mean(GhStd(:))/mean(BGStd(:)))}, 'Parent',ax)
	xlabel(ax, 'Image intensity (a.u.)')
	ylabel(ax, 'Nr of voxels')
	title(ax, Title)
	
	ax = subplot(4,1,3, 'Parent',HG);
	BGAvg = BGAvg';
	plot(ax, BGAvg(:))
	xlim(ax, [1 numel(BGAvg)])
	xlabel(ax, 'Slices(t)')
	ylabel(ax, 'mean(Noise)')
	hold on
	plot(ax, [1 numel(BGAvg)], Thresh * [1 1], '--')
	if Spikes(1)>0
		text(0.75*numel(BGAvg), max(BGAvg(:)), sprintf('Nr of spikes: %d', Spikes(1)), 'Parent',ax)
	end
	
	ax = subplot(4,1,4, 'Parent',HG);
	% BGStd = BGStd';
	plot(ax, linspace(1,size(BGStd,2),numel(BGStd)), BGStd(:))
	xlim(ax, [1-eps size(BGStd,2)])							% Robust against single image input
	set(ax, 'XTick', unique(round(get(ax,'XTick'))))
	xlabel(ax, 'Slices(z)')
	ylabel(ax, 'std(Noise)')
	
end
