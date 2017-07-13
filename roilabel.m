function CMask = roilabel(ROIMask_)

% ROILABEL finds ROI clusters in a Mask and labels them with increasing order
% (1 for ROI1, 2 for ROI2, etc). Typically these ROIs are white
% matter lesions that are segmented using MRICro.
%
% FUNCTION CMask = ROILABEL(ROIMask)
%
% INPUT
%   ROIMask	- String array of ROI-file or nD binary Mask volume
%
% OUPUT
%   LMask	- Labeled ROI mask
%
% See also: BWLABEL
%
% -- Marcel, 04-02-2011.

if nargin<1 || isempty(ROIMask_)
	ROIMask_ = spm_select(1, 'image', 'Select unlabeled Mask file');
end
if ischar(ROIMask_)
	disp('-------------------------------------------------------')
    disp(['Reading ' ROIMask_])
    ROIMask_ = round(spm_read_vols(spm_vol(ROIMask_)));	% => partial voluming / precision
end

% Add zero-sheet to prevent errors in region growing
MaskSz  = [1 1 1];
MaskSz(1:ndims(ROIMask_)) = size(ROIMask_);
ROIMask = zeros(MaskSz + [2 2 2]);
ROIMask(2:end-1, 2:end-1, 2:end-1) = ROIMask_;

% Start identifying clusters and give them an increasing nr
%
% NB: would have been easier and up to 4x faster when using:
% ROIMask = bwlabeln(ROIMask);
%
% My coding effort below was a nice programming excercise though :-)
%
if size(ROIMask,3)==1
	[X Y Z] = meshgrid(-1:1, -1:1, 0);			% 2D neighbouring voxels
else
	[X Y Z] = meshgrid(-1:1, -1:1, -1:1);		% 3D neighbouring voxels
end
ClusNr = 1;
NZVox  = find(ROIMask);
while ~isempty(NZVox)
	Seeds = 1;
	Clus  = NZVox(Seeds);						% Take the first non-zero voxel
	ROIMask(Clus) = ClusNr+1;					% Update the value of this voxel
	while Seeds
		[I J K] = ind2sub(size(ROIMask), Clus(Seeds));
		OldClus = Clus;
		for i = 1:length(X(:))
			% Region growing
			ClusN = sub2ind(size(ROIMask), I+X(i), J+Y(i), K+Z(i));
			ClusN = ClusN(ROIMask(ClusN)==1);	% Find non-zero neighbouring voxels
			Clus  = [Clus; ClusN];				% Update cluster-list
			ROIMask(ClusN) = ClusNr+1;			% Update non-zero neighbouring voxels
		end
		Seeds = [length(OldClus)+1:length(Clus)];
	end
	ClusNr = ClusNr + 1;
	% Keep the remaining non-zero voxels
	NZVox = setdiff(NZVox, Clus);
end

% Subtract the original ROIMask to get numbering right
ROIMask(ROIMask>0) = ROIMask(ROIMask>0) - 1;
CMask = ROIMask(2:end-1, 2:end-1, 2:end-1);