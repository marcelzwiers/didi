function varargout = dd_initcnode(Job)

persistent InitCfg

if nargin<1 || Job.ParallelBox.Val
	if isempty(InitCfg)
		TmpDir = tempname;			% Avoid tempdir as it may slow matlab down a lot and suffer from (spm_print) concurrency problems
		if ~exist(TmpDir, 'dir')
			mkdir(TmpDir)
		end
		cd(TmpDir)					% Dump spm_print output here (e.g. from spm_coreg)
		spm_jobman('initcfg')		% Seems necessary (p2p does not transfer this)
		[S LWarn] = mywarning('Off', 'SPM:noDisplay');
		spm_figure('GetWin', 'Graphics');
		spm_figure('GetWin', 'Interactive');
		mywarning(S, LWarn)
		spm_get_defaults('cmdline', true)
		InitCfg = true;
	end
end

if nargout
	varargout{1} = spm_figure('FindWin', 'Graphics');
	varargout{2} = spm_figure('FindWin', 'Interactive');
end
