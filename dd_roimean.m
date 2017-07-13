function varargout = dd_roimean(Inds, ROIs, LabelNrs)

% FUNCTION VARARGOUT = DD_ROIMEAN(IndicesFiles, ROIMasks, LabelNrs)
%
% DD_ROIMEAN computes the (weighted) mean of diffusion tensor derivatives of
% interest (e.g. FA or MD) over specified regions-of-interest.
%
% INPUT
%	IndicesFiles - Character or cell array of file-names of ascii-files (e.g.
%				   dd_2011Jul27_Job_fa.txt) that list the full pathnames of
%				   individual diffusion tensor derivative volumes of interest
%				   (e.g. of the fa-volumes). This argument can be left empty for
%				   interactive use.
%	ROIMasks	 - Character or cell array of file-names of the masking /
%				   weighting volume(s) that are multiplied with the tensor
%				   derivative volumes. These volumes must be in the same space
%				   as the tensor derivative volumes. This argument can be left
%				   empty for interactive use.
%	LabelNrs	 - Array of label numbers that constitute your ROI in ROIMasks
%				   (such as when using JHU ICBM-DTI-81 White-Matter Labels)
%
% OUTPUT
%	VARARGOUT	 - An output variable is given for each IndicesFile. The output
%				   variables are nr-rows-of-indicesfile x nr-of-roimasks matrices
%
% EXAMPLES
%	>>dd_roimean
%	>>[FA MD] = dd_roimean;
%	>>[FA MD MO] = dd_roimean(['../Jobs/dd_2011Jul27_Job_fa.txt'; ...
%		'../Jobs/dd_2011Jul27_Job_ad.txt';'../Jobs/dd_2011Jul27_Job_mo.txt'], ...
%		{'../Masks/MyAmygd.nii', fullfile(spm('dir'),'toolbox','Anatomy','PMaps','Fiber_forn.img')});
%
% Marcel, 28-07-2011

if nargin<1 || isempty(Inds)
	if nargout
		NFiles = nargout;
	else
		NFiles = Inf;
	end	
	Inds = spm_select(NFiles, 'mat', 'Select tensor indices filelist-file(s) (e.g. *Job_fa.txt etc)');
end
if nargin<2 || isempty(ROIs)
	ROIs = spm_select(Inf, 'image', 'Select ROI-mask(s)');
end
if nargin<3
	LabelNrs = [];
end

ROIs	= spm_vol(ROIs);
ROIVols = spm_read_vols(ROIs);
if ~isempty(LabelNrs)
	ROIVols0 = zeros(size(ROIVols));
	for n = 1:numel(LabelNrs)
		ROIVols0(ROIVols==LabelNrs(n)) = 1;
	end
	ROIVols = ROIVols0;
end
ROIVols(isnan(ROIVols)) = 0;			% Account for NaNs => sum()
for n = 1:size(ROIVols,4)				% Normalize the ROI-masks => weights
	if nargin<2
		fprintf('Column{%d}:\t%s\n', n, ROIs(n).fname)
	end
	ROIVols(:,:,:,n) = ROIVols(:,:,:,n) / sum(sum(sum(ROIVols(:,:,:,n))));
end

for n = 1:size(Inds,1)
	if nargin<1
		fprintf('\nVARARGOUT{%d}:\t%s\n', n, Inds(n,:))
	end
	FID     = fopen(strtrim(Inds(n,:)));
	IndVols = textscan(FID, '%s');
	fclose(FID);
	IndVols = spm_read_vols(spm_vol(char(IndVols{1})));
	IndVols(isnan(IndVols)) = 0;		% Account for NaNs => sum()
	varargout{n} = NaN(size(IndVols,4), size(ROIVols,4));
	for m = 1:size(ROIVols,4)			% Compute the weighted sum (= weighted mean)
		varargout{n}(:,m) = sum(sum(sum(bsxfun(@times, IndVols, ROIVols(:,:,:,m)))));
	end
	if ~nargout
		disp(varargout{n})
	end
end

if ~nargout, varargout = {}; end