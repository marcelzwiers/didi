function [NrVox V] = dd_eigvalcorr(DImgs, Mask, PreFix, CVox)

% FUNCTION [NrVox cDImgs] = dd_eigvalcorr(DImgs, Mask, PreFix, CVox)
%
% Correct the diffusion tensor bij setting the (non-physical) negative
% eigenvalues to zero.
%
% INPUT
%   DImgs  - String or cell array containg DT filenames
%            e.g. {Dxx_* Dxy_* Dxz_* Dyy_* Dyz_* Dzz_*}
%	Mask   - Filename of Mask or Boolean volume to limit NrVox (NB: This does
%			 not effect cDImgs)
%   PreFix - The output is written to pre-fixed files
%            (leave empty to overwrite).
%   CVox   - Saves a volume with ones at the position of the corrected
%            voxels if true (default is false)
%
% OUTPUT
%	NrVox  - Number of (in-mask) voxels for which the eigenvalues were corrected
%   cDImgs - Header information of corrected DT filenames
%
% Marcel, 10-01-2008.

if nargin<1 || isempty(DImgs)
    DImgs = spm_select(6, '((dt[1-6])|(D[x-z][x-z]))_.*(\.img$|\.nii$)');
end
DImgs = spm_vol(char(DImgs));
if numel(DImgs) ~= 6
    error('Wrong number of tensor elements');
end
if nargin<2 || isempty(Mask)
	Mask = true(DImgs(1).dim);
elseif ischar(Mask) || iscellstr(Mask)
	Mask = logical(spm_read_vols(spm_vol(char(Mask))));
end
if nargin<3, PreFix = ''; end
if nargin<4, CVox = false; end

% Load the diffusion tensor volumes
fprintf('\nSetting negative eigenvalues to zero: %s\n', DImgs(1).fname)
DVol = permute(spm_read_vols(DImgs), [4 1 2 3]);	% Using first index avoids 1 for loop
if CVox
    DSz  = size(DVol);
    CVol = false(DSz(2:4));
end

% Loop over slices/voxels to correct negative eigenvalues
Count = 0;
fprintf('Slice: %03d/%03d', 0, size(DVol,2));
spm_progress_bar('Init', size(DVol,2), 'Setting negative eigenvalues to zero', 'Slices completed');
for x = 1:size(DVol,2)					% Loop over all y-z planes
    fprintf([repmat('\b',[1 7]) '%03d/%03d'], x, size(DVol,2));
	spm_progress_bar('Set', x)
    for n = 1:numel(DVol(1,x,:))		% Loop over all voxels in y-z plane
        if all(isfinite(DVol(:,x,n)))	% Input to eig must not contain NaN or Inf
            [Vec Val] = eig([DVol(1,x,n) DVol(2,x,n) DVol(3,x,n)
                             DVol(2,x,n) DVol(4,x,n) DVol(5,x,n)
                             DVol(3,x,n) DVol(5,x,n) DVol(6,x,n)]);
            Lbd = diag(Val);
            if any(Lbd<=0)              % Val = matrix of eigenvalues (see eig)
                Lbd(Lbd<=0) = eps*randi(100,sum(Lbd<=0),1); % Use Val<=0 to avoid division-by-zero warnings (e.g. in FA), use randi to avoid spurious FA = 0 values causing troubles in TBSS_1_preproc
                D = Vec*diag(Lbd)/Vec;	% Vec = matrix of eigenvectors (see eig)
                DVol(:,x,n) = [D(1,1); D(1,2); D(1,3); D(2,2); D(2,3); D(3,3)];
				if Mask(x,n)			% Only count the in-mask corrections
					Count = Count + 1;
				end
                if CVox
                    CVol(x,n) = true;
                end
            end
        end
    end
end
fprintf('\nNr of (in-mask) corrected voxels: %d/%d\n', Count, sum(Mask(:)))

% Save the diffusion tensor volumes
for n = 1:6
    V(n) = spm_write_vol(prepend(DImgs(n), PreFix), shiftdim(DVol(n,:,:,:)));
end

% Save the correction volume
if CVox
    CHdr         = prepend(DImgs(1), 'CVox_');
    CHdr.descrip = 'Corrected negative eigenvalues';
    CHdr.dt(1)   = spm_type('uint8');   % Set datatype to uint8
    spm_write_vol(CHdr, CVol);
end

if nargout
	NrVox = [Count sum(Mask(:))];
end


%% ---------------------- END ----------------------

function Vol = prepend(Vol, PreFix)

for n = 1:numel(Vol)
    [Pth, Nm, Xt] = fileparts(Vol(n).fname);
    Vol(n).fname  = fullfile(Pth, [PreFix Nm Xt]);
end
Vol = rmfield(Vol, 'pinfo');
