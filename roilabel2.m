function [CMask ClusNr] = roilabel2(ROIMask_, Conn)

% ROILABEL2 finds ROI clusters in a Mask and labels them with increasing order
% (1 for ROI1, 2 for ROI2, etc). Typically these ROIs are 2D DTI artefacts.
%
% FUNCTION [CMask NUM] = ROILABEL(ROIMask, N)
%
% INPUT
%   ROIMask	- 2D binary Mask volume
%	N		- Connectivity (4 or 8)
%
% OUPUT
%   LMask	- Relabeled Mask volume
%	NUM		- The number of connected object found in ROIMask
%
% See also: ROILABEL, ROILABELN, BWLABEL, DD_PATCH.
%
% -- Marcel, 04-02-2011.

% Add zero-sheet to prevent errors in region growing
MaskSz_ = size(ROIMask_);
ROIMask = zeros(MaskSz_ + [2 2]);
ROIMask(2:end-1, 2:end-1) = ROIMask_;
MaskSz  = size(ROIMask);

% Start identifying clusters and give them an increasing nr
%
% NB: would have been easier and ~10-50x faster when using:
% ROIMask = bwlabel(ROIMask);
%
% My coding effort below was a nice programming excercise though :-)
%
[X Y] = meshgrid(-1:1, -1:1);					% 2D neighbouring voxels
if nargin<2 || Conn==4;
	X([1 3 7 9]) = [];
	Y([1 3 7 9]) = [];
end
ClusNr = 0;
NZVox  = find(ROIMask);
while ~isempty(NZVox)
	ClusNr = ClusNr + 1;
	Seeds  = 1;
	Clus   = NZVox(Seeds);						% Take the first non-zero voxel
	ROIMask(Clus) = ClusNr+1;					% Update the value of this voxel
	while Seeds
		[I J] = ind2sub(MaskSz, Clus(Seeds));
		OldClus = Clus;
		for i = 1:length(X(:))
			% Region growing
			ClusN = mysub2ind(MaskSz, I+X(i), J+Y(i));
			ClusN = ClusN(ROIMask(ClusN)==1);	% Find non-zero neighbouring voxels
			Clus  = [Clus; ClusN];				% Update cluster-list
			ROIMask(ClusN) = ClusNr+1;			% Update non-zero neighbouring voxels
		end
		Seeds = [length(OldClus)+1:length(Clus)];
	end
	% Keep the remaining non-zero voxels
	NZVox = setdiff(NZVox, Clus);
end

% Subtract the original ROIMask to get numbering right
ROIMask(ROIMask>0) = ROIMask(ROIMask>0) - 1;
CMask = ROIMask(2:end-1, 2:end-1);


function ndx = mysub2ind(siz,varargin)

% Stripped down for speed

%Compute linear indices
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
for i = 1:length(siz)
    ndx = ndx + (varargin{i}-1)*k(i);
end