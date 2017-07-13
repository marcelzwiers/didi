function Outliers = dd_outliers(RefImage, Images, Range)

% DD_OUTLIERS is an interactive tool to detect images that are outside a given
% range for a given location.
%
% INPUT
%	RefImage - 
%	Images	 -
%	Range	 -
%
% OUTPUT
%	Outliers -
%
% Marcel, 8-9-2016

% Defaults
if nargin<1 || isempty(RefImage)
	RefImage = spm_select(1, 'image', 'Select the reference image');
end
if nargin<2 || isempty(Images)
	Images = spm_select(Inf, 'image', 'Select the image files for testing');
end
Images = cellstr(Images);
if nargin<3 || isempty(Range)
	Range = [-1 1];
end
TmpDir = tempname;
for n = 1:numel(Images)
	[~, Nme Ext] = fileparts(Images{n});
	if strcmp(Ext, '.gz')
		fprintf('Unzipping to: %s\n', fullfile(TmpDir,Nme))
		gunzip(Images{n}, TmpDir)
		Images{n} = fullfile(TmpDir,Nme);
	end
end

% Get the images
Hdrs   = spm_vol(char(Images));
RefHdr = spm_vol(RefImage);

% Display the reference image
spm_check_registration(RefHdr)
spm_orthviews('Interp', 0)							% No interpolation = better visibility of artefacts
spm_orthviews('context_menu','orientation',2)		% Use voxel space (1st image)

% Let the user choose a location and get the outliers
Finish = 
while ~Finish
	Pos		 = RefHdr.mat \ [spm_orthviews('Pos'); 1];
	Y		 = spm_get_data(Hdrs, round(Pos(1:3)));
	Outliers = Hdrs(Y<=Range(1) | Y>=Range(2));
end

% Clean up the temporary unzipped files
if exist(TmpDir,'dir')
	rmdir(TmpDir,'s')
end
