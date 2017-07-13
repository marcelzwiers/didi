function varargout = dd_wregress(rFNames, rWTot, PreFix)

% FUNCTION DVol = dd_wregress(rFNames, rWTot, PreFix)
%
% INPUT:
%	rFNames	- The filenames of the resliced DW images
%	rWTot	- The resliced tzyx-weights (n.b. not directly from PATCH.mat)
%	PreFix	- Prefix for the saved D##-files (if nargout==0)
%
% OUTPUT
%	DVol	- 4D matrix with the D## data (the same dimflip as in dd_dt_regress; see dd_rotategradients)
%
% Marcel, 6-8-2010

% This permutation order is convenient for the regression calculations
tzyx = [4 3 2 1];
B0	 = 50;											% Threshold for selection of b0 images

% Parse the input variables
if nargin<3 || isempty(PreFix)
	PreFix = {'Dxx_PATCH_', 'Dxy_PATCH_', 'Dxz_PATCH_', 'Dyy_PATCH_', 'Dyz_PATCH_', 'Dzz_PATCH_'};
end
rFNames = char(rFNames);
if ischar(rWTot);
	rWTot = permute(single(spm_read_vols(spm_vol(rWTot))), tzyx);
end

% Get some DTI info and data
disp('Loading DW volumes & weights...')
spm('FigName', 'PATCH estimation');
for n = 1:size(rFNames,1)
	D(n) = dti_get_dtidata(rFNames(n,:));
end
DWISel = [D(:).b] > B0;								% The DWI files
q      = vertcat(D(DWISel).g);						% Matrix containing DTI gradient directions in world-space
b      = vertcat(D.b);								% Collumn vector containing all b-values
VolHdr = spm_vol(rFNames);
Vol	   = permute(single(spm_read_vols(VolHdr)), tzyx);
VolSz  = size(Vol);
Vol(Vol<1) = 1;										% Avoid log(0) = -Inf and log(-#) = complex

%-- Compute the WLS tensor coefficients
% S   = S0*exp(-BD) + E
% [y  = log(S0/S)]
% [D  = [Dxx Dyy Dzz Dxy Dxz Dyz]]
% BD  = yy = y-e
% WBD = Wy - We
% D   = inv((WB)'WB)(WB)'Wy   (NB: the expectation value of x*e is zero)

disp('Estimating the tensor coefficients using PATCH...')
B = repmat(b(DWISel), [1 6]) .* [q.^2 2*q(:,1).*q(:,2) 2*q(:,1).*q(:,3) 2*q(:,2).*q(:,3)];
B = B(:, [1 4 5 2 6 3]);					% Make D alphabetical: [Dxx Dxy Dxz Dyy Dyz Dzz] 
y = log(bsxfun(@rdivide, mean(Vol(~DWISel,:),1), Vol(DWISel,:)));

% Compute the diffusion tensor per voxel
NrUpd  = round(prod(VolSz(2:end))/250);		% Divide display waitbar in +/-250 pieces
spm_progress_bar('Init', prod(VolSz(2:end)), 'Weighted regression', 'Voxels completed');
DVol   = zeros([6 VolSz(2:end)]);
CondNr = 10*eps(class(y));
for n = 1:prod(VolSz(2:end))
	W  = rWTot(:,n);
	WB = bsxfun(@times, W, B);				% bsxfun is faster than repmat
	if rcond(WB'*WB) < CondNr				% Test if the matrix inverse can be properly computed (i.e. if we can make a WLS estimate)
		WB = B;								% If not then revert to OLS
	end
	DVol(:,n) = (WB'*WB)\WB'*(W.*y(:,n));
	if rem(n, NrUpd)==0
		spm_progress_bar('Set', n)
	end
end

%-- Save the results to the space of the first DTI-volume
Hdr		  = rmfield(VolHdr(1), 'pinfo');	% TODO: check if this is a b0-volume?
Hdr.dt(1) = spm_type('float32');
[FPath FName Ext] = fileparts(Hdr.fname);
disp('Saving estimated tensor coefficients:')
for n = 1:size(DVol,1)
	Hdr.fname = fullfile(FPath, [PreFix{n} FName Ext]);
	disp([PreFix{n} FName Ext])
	spm_write_vol(Hdr, ipermute(DVol(n,:,:,:), tzyx));
end
if nargout
	varargout{1} = ipermute(DVol, tzyx);
end


%% -- END --

% % function Vol = bireslice(Vol, Mi)
% % 
% % VolSz = size(Vol);
% % spm_progress_bar('Init', VolSz(1), 'Reslicing', 'Volumes');
% % for ti = 1:VolSz(1)
% % 	spm_progress_bar('Set', ti)
% % 	Volxyz = permute(Vol(ti,:,:,:), [4 3 2 1]);			% Vol = Vol(t,z,y,x)
% % 	for zi = 1:VolSz(2)
% % 		ni = min(size(Mi,4), zi);						% ni is either 1 or zi (Andersson)
% % 		if any(any(abs(Mi(:,:,ti,ni)-eye(4)) > 1e-12))	% Account for rounding errors
% % 			% Reslice the xyz-volume using sinc interpolation; reslicing to realigned space => inverse rotation
% % 			Vol(ti,zi,:,:) = ipermute(spm_slice_vol(Volxyz, Mi(:,:,ti,ni)\spm_matrix([0 0 zi]), ...
% % 													VolSz([4 3]), -1), [2 1]);
% % 		end
% % 	end
% % end
