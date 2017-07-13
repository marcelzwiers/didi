function [D DimFlip] = dd_rotategradients(D, M, DimFlip)

% FUNCTION [D DimFlip] = dd_rotategradients(D, M, DimFlip)
%
% This function exists because Volkmar used a LAS instead of the RAS Tailarach
% coordinate system, i.e. he did not flip the sign of the x-component of the
% gradient vectors; see
% http://sourceforge.net/projects/spmtools/forums/forum/485380/topic/2125699 and
% spm_dicom_convert.m and dicominfodti.m
%
% INPUT:
%	D		- Structure from dti_get_dtidata (with D.g field) or nx3 gradient matrix
%	M		- Rotation that needs to be applied to the gradient vectors
%	DimFlip	- The x-, y- and z-fip that was used to bring the vectors into world coordinates
%
% OUTPUT:
%	D		- Structure/matrix with rotated gradient vectors
%	DimFlip	- The x-, y- and z-fip that was used to bring the vectors into world coordinates
%
% Marcel, 25-08-2010.

if nargin<3 || isempty(DimFlip)
	DimFlip = [-1 1 1];
end
if nargin<2 || all(all(M==eye(4)))	% We're not gonna do rotation, just return the DimFlip
	if nargin<1
		D = [];
	end
	return
end

Par = spm_imatrix(M);
R	= spm_matrix([0 0 0 Par(4:6)]);		% Remove translation & scaling (if any) for pure rotation
if isstruct(D)
	for n = 1:numel(D)
		D(n).g = DimFlip .* (R(1:3,1:3) * (DimFlip .* D(n).g)')'; % Do x-flip, rotate and x-flip back
	end
else
	D = (diag(DimFlip) * R(1:3,1:3) * diag(DimFlip) * D')';
end