function Job = dd_basicproc_cleanup(Job, SubjNrs, SeriesNrs, Mode)

% DD_BASICPROC_CLEANUP is a key internal function of dd_basicproc that
% is made available externally to allow distributed computing.
%
% FUNCTION Job = DD_BASICPROC_CLEANUP(Job, SubjNr, SeriesNr, Mode)
%
% Marcel, 12-7-2011
%
% See also: DD_BASICPROC

if ~Job.CleanUpBox.Val
	return
end
if nargin<2 || isempty(SubjNrs)
	SubjNrs = 1:size(Job.Nifti,1);
end
if nargin<3 || isempty(SeriesNrs)
	SeriesNrs = 1:size(Job.Nifti,2);
end
if nargin<4 || isempty(Mode)
	Mode = {'Denoise', 'PATCH', 'DTENSOR'};
end
disp('Cleaning up intermediate denoised-, PATCH- or DW-images...')

for SubjNr = SubjNrs(:)'
	for SeriesNr = SeriesNrs(:)'
		try
			
		% Delete all denoised files
		if any(strcmp('Denoise', Mode)) && ~strcmp(Job.DenoisingMenu.Str{Job.DenoisingMenu.Val}, 'none')
			delete(fullfile(Job.Nifti(SubjNr,SeriesNr).Path, [Job.DenoisingMenu.Str{Job.DenoisingMenu.Val} '_*']))
		end

		% Delete all PATCH files
		if any(strcmp('PATCH', Mode))
			delete(fullfile(Job.Nifti(SubjNr,SeriesNr).Path, 'PATCH*'))
		end	

		% Delete everything except the original Nifti- and DICOM-files + all tensor-files + export-files
		% NB: Be careful not to delete original Nifti-files (PreFix='')
		if any(strcmp('DTENSOR', Mode)) && ...
		   ~(~(Job.T1Box.Val && Job.PEUnwarpBox.Val) && strcmp(Job.RealignMenu.Str{Job.RealignMenu.Val}, 'none'))
			IntImgs = regexprep([cellstr(Job.Output(SubjNr,SeriesNr).WImgs)' Job.Output(SubjNr,SeriesNr).TgtImgs], '\.img$|\.nii$', '.*');
			delete(IntImgs{~cellfun(@isempty,IntImgs)}, fullfile(Job.Nifti(SubjNr,SeriesNr).Path, 'artc_*'), fullfile(Job.Nifti(SubjNr,SeriesNr).Path, 'meanpw*'))
			Job.Output(SubjNr,SeriesNr).TgtImgs = {};
			Job.Output(SubjNr,SeriesNr).WImgs   = '';
		end
		
		fprintf('.')
		
		end
	end
end
fprintf('\n')
