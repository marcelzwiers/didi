function Job = dd_bids(BIDSDir, Sub, Ses, DidiDir)

% Job = DD_BIDS(BIDSDir, Sub, Ses, DidiDir)
%
% Creates a shadow directory in bids/derivatives/didi/.. with unzipped copies of the DWI data

fprintf('Building the BIDS-map\n')
if nargin<1 || isempty(BIDSDir)
	BIDSDir = spm_select(1,'dir','Select the BIDS-directory');
elseif BIDSDir(1) ~= filesep					% Turn a realtive path into an absolute one. TODO: Make this OS-portable
	BIDSDir = fullfile(pwd, BIDSDir);
end
BIDS = spm_BIDS(BIDSDir)
if nargin<2 || isempty(Sub)
	Sub	= spm_file(spm_BIDS(BIDS, 'subjects', 'type','dwi'), 'prefix','sub-');
else
	Sub	= spm_file(spm_BIDS(BIDS, 'subjects', 'sub',Sub), 'prefix','sub-');
end
if nargin<3 || isempty(Ses)
	Ses	= spm_file(spm_BIDS(BIDS, 'sessions', 'type','dwi'), 'prefix','ses-');
else
	Pre       = strncmp('ses-',Ses,4);
	Ses(~Pre) = spm_file(Ses(~Pre), 'prefix','ses-');
end
if isempty(Ses)
    Ses = {''};
end
if nargin<4 || isempty(DidiDir)
	DidiDir = fullfile(BIDSDir, 'derivatives','didi');
elseif DidiDir(1) ~= filesep					% Turn a realtive path into an absolute one. TODO: Make this OS-portable
	DidiDir = fullfile(pwd, DidiDir);
end

fprintf('Copying the data\n')
spm_mkdir(DidiDir, Sub, Ses, 'dwi');
spm_mkdir(DidiDir, Sub, Ses, 'anat');
[DWI_dst, T1_dst] = deal(cell(numel(Sub), numel(Ses)));
for n = 1:numel(Sub)
	for m = 1:numel(Ses)
		fprintf('-> %s\n', fullfile(DidiDir, Sub{n}, Ses{m}))
		DWI_src		 = spm_BIDS(BIDS, 'data', 'sub',Sub{n}, 'ses',Ses{m}, 'type','dwi');
		T1_src		 = spm_BIDS(BIDS, 'data', 'sub',Sub{n}, 'ses',Ses{m}, 'type','T1w');
		assert(size(T1_src,1) < 2, 'More than 1 T1w image found, do not know which one to choose...')
		DWI_dst{n,m} = spm_file(DWI_src, 'path', fullfile(DidiDir,Sub{n},Ses{m},'dwi'));
		T1_dst{n,m}  = spm_file(T1_src,  'path', fullfile(DidiDir,Sub{n},Ses{m},'anat'));
		spm_copy(DWI_src, spm_file(DWI_dst{n,m},'path'), 'nifti',true, 'gunzip',true)
		spm_copy(T1_src,  spm_file(T1_dst{n,m}, 'path'), 'nifti',true, 'gunzip',true)
		% spm_copy does not (yet?) support bval/bvec file-extensions (neither does it consider a .nii.gz extension as nifti)
		spm_copy(spm_file(spm_file(DWI_src,'ext',''),'ext','.bv*'), spm_file(DWI_dst{n,m},'path'))
	end
end

% Create the Job-structure
fprintf('\nBuilding the Job structure ')
for n = 1:numel(Sub)
	for m = 1:numel(Ses)
		fprintf('.')
		Job.Nifti(n,m).Path	 = spm_file(DWI_dst{n,m}{end},'path');
		Job.Nifti(n,m).Files = spm_file(spm_file(DWI_dst{n,m},'basename'), 'ext','.nii');
		if numel(Job.Nifti(n,m).Files) ~= 1
			disp(char('Found:', Job.Nifti(n,m).Files{:}))
			for File = Job.Nifti(n,m).Files(1:end-1)'
				delete(char(spm_file(File, 'path',Job.Nifti(n,m).Path, 'ext','.*')))
			end
			Job.Nifti(n,m).Files = Job.Nifti(n,m).Files(end);		% It's most likely that the last one is the best one?
			warning('Cannot (yet) deal with multiple DWI runs, will use: %s', char(Job.Nifti(n,m).Files))
		end
		Job.Nifti(n,m).T1Text = char(spm_file(spm_file(T1_dst{n,m}, 'ext',''),'ext','.nii'));
		if size(Job.Nifti(n,m).T1Text,1) > 1
			disp(char('Found:', Job.Nifti(n,m).T1Text))
			Job.Nifti(n,m).T1Text = deblank(Job.Nifti(n,m).T1Text(end,:));	% It's most likely that the last one is the best one?
			warning('There should exists one (best) T1-image, will use: %s', Job.Nifti(n,m).T1Text)
		end
	end
end
Job.Ver			 = 'dd_basicproc 4.0';		% This should always be updated and adapted to the latest version
Job.DICOM		 = [];						% Forget about the DICOMs
Job.DataMenu.Str = {'DICOM';'Nifti'};
Job.DataMenu.Val = 2;						% We have nifti-files
Job.SubjMenu.Str = cellfun(@num2str, num2cell(1:numel(Job.Nifti))', 'UniformOutput',false);
Job.SubjMenu.Val = 1;						% Irrelevant
Job.T1Box.Val	 = ~isempty(T1_src);		% We have a T1 image and want to be in coregister with it
Job.BETMenu.Str	 = '';						% Needed

% Save the Job to the derivative-folder (dd_basicproc will store more output here) and launch the GUI
JobName = fullfile(DidiDir, ['dd_' datestr(now,'yyyymmmdd') '_Job']);
fprintf('\nSaving Job-info to: %s\n', JobName)
save(JobName, 'Job')
H = dd_basicproc;
dd_basicproc('LoadButton_Callback', H, [], guidata(H), Job)
