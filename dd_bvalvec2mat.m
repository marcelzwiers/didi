function D = dd_bvalvec2mat(DWIFiles, BValFile, BVecFile, Coord)

% FUNCTION D = dd_bvalvec2mat(DWIFiles, BValFile, BVecFile, Coord)
%
% DD_BVALVEC2MAT reads b-values and b-vectors from disk and stores them in the
% .mat (LAS) file format of the SPM diffusion toolbox of Volkmar Glauche.
%
% INPUT
%	DWIFiles - Character or cell array with the names of the DWI files. The order
%			   must match the order in BValFile and BVecFile
%	BValFile - Filename of the ascii-file containing the b-values
%	BVecFile - Filename of the ascii-file containing the b-vectors. The format
%			   must be [n x 3] or [3 x n].
%	Coord	 - The b-vectors can either be represented in RAS world coordinates
%			   (Coord = 'World') or in voxel coordinates (Coord = 'Voxel'). If
%			   empty then Coord = 'Voxel' (as in FSL).
%
% OUTPUT
%	D		 - A structure containing the stored b-vec and b-val information
%
% See also: DICOMINFODTI, DD_BASICPROC_EXPORTDWI, DD_ROTATEGRADIENTS, DTI_GET_DTIDATA
%
% Marcel, 6-3-2015

% Defaults
if nargin<1 || isempty(DWIFiles)
	DWIFiles = spm_select(Inf,'image','Select the DWI files');
else
	DWIFiles = char(DWIFiles);
end
if nargin<2 || isempty(BValFile)
	BValFile = spm_select('FPList', fileparts(DWIFiles(1,:)), 'bval.txt|.*bval$');	% Make an educated guess
	if isempty(BValFile)
		BValFile = spm_select(1,'any','Select the file containg the b-value information',{},fileparts(DWIFiles(1,:)));
	end
end
if nargin<3 || isempty(BVecFile)
	BVecFile = strrep(BValFile, 'bval', 'bvec');									% Make an educated guess
	if ~exist(BVecFile,'file')
		BVecFile = spm_select(1,'any','Select the file containg the b-vector information',{},fileparts(BVecFile));
	end
end
if nargin<4 || isempty(Coord)
	warning('No coordinate system specified, assuming voxel coordinates (as in FSL)')
	Coord = 'Voxel';	% FSL & Camino
end

% Get the bval- and bvec information
BVec = load(BVecFile);
Sz	 = size(BVec);
if ~any(Sz==3) || numel(Sz)~=2
	error('The b-vector file (%s) has wrong size [%d x %d] instead of [n x 3] or [3 x n]', BVecFile, Sz)
elseif Sz(1)==3
	BVec = BVec';
end
BVal = load(BValFile);
if numel(BVal)~=size(BVec,1)
	error('The b-value file (%s) has wrong size [%d] instead of [%d]', BValFile, numel(BVal), size(BVec,1))
end
if size(DWIFiles,1) ~= numel(BVal)
	error('The number of DWI-files (%d) does not match with the number of bval/bvec elements (%d; See e.g. %s)', size(DWIFiles,1), numel(BVal), BValFile)
end

% Map the vectors from RAS to LAS coordinates, i.e. include Volkmar's x-flip error! (see dd_rotategradients)
if strcmpi(Coord, 'voxel')							% Transform the DW directions (vectors) from voxel- to world RAS coordinates
	for n = 1:size(DWIFiles,1)
		RPar	  = spm_imatrix(spm_get_space(DWIFiles(n,:)));
		R		  = spm_matrix([0 0 0 RPar(4:6)]);	% Remove translation & scaling (if any) for pure rotation
		BVec(n,:) = BVec(n,:) * R(1:3,1:3)';		% = (R(1:3,1:3)*BVec(n,:)')'; See also: dd_basicproc_exportdwi
	end
end
[D Flip] = dd_rotategradients;
BVec	 = bsxfun(@times, BVec, Flip);				% RAS -> LAS

% Store the DWI info to the mat file
for n = 1:size(DWIFiles,1)
	D(n).g	 = BVec(n,:);
	D(n).b	 = BVal(n);
	D(n).mat = spm_get_space(DWIFiles(n,:));
	dti_get_dtidata(DWIFiles(n,:), D(n));
end
