function dd_disp_grad(Files, Options)

% Visualise gradient direction vectors
%
% A wrapper around DTI_DISP_GRAD to visualise a set of gradient direction
% vectors. This is useful to compare gradient directions between measurements
% which may differ due to changes in measurement protocols, motion correction
% etc.
%
% See also: DTI_DISP_GRAD
%
% Marcel, 28-6-2015

if ~exist('dti_disp_grad','file')
	addpath('/home/common/matlab/spm12', '/home/common/matlab/spm12/toolbox/Diffusion/Helpers', '/home/common/matlab/spm12/toolbox/Diffusion/Visualisation')
end
if nargin<1 || isempty(Files)
	Files = spm_select(Inf, 'image');
end
if nargin<2 || isempty(Options)
	Options = 'vd';
end

Bch.srcimgs = cellstr(Files);
Bch.option	= Options;
Bch.fgbgcol = [0 0.2 0.2];
Bch.axbgcol = [0 0 0];
Bch.res		= 50;

dti_disp_grad(Bch)
close(spm_figure('FindWin','Interactive'))
set(gca, 'CLim', max(get(gca,'CLim')) * [0.5 1])
