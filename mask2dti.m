function mask2dti(DefFile, Masks, OutputDir)

if nargin<3 || isempty(OutputDir)
	OutputDir = spm_select(1, 'dir', 'Select the output directory');
end
if nargin<2 || isempty(Masks)
	Masks = spm_select(Inf, 'image', 'Select the mask file(s)');
end
if nargin<1 || isempty(DefFile)
	DefFile = spm_select(1, 'image', 'Select the inverse deformation file', [], [], '^iy_.*');
end

spm_jobman('initcfg')

[Tag Job] = spm_jobman('harvest', spm_cfg_defs);
Job.comp{1}.def{1} = DefFile;
Job.interp		   = 1;

Job.savedir = struct('saveusr',{{OutputDir}});
for n = 1:size(Masks,1)
	Job.fnames{n,1} = deblank(Masks(n,:));
end
Warp{1}.spm.util.(Tag) = Job;
spm_jobman('run', Warp);
