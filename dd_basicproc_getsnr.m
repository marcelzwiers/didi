function OK = dd_basicproc_getsnr(Job, SubjNr, LogName)

% DD_BASICPROC_GETSNR is a key internal function of dd_basicproc that is made
% available externally to allow distributed computing.
%
% Marcel, 16-4-2014
%
% See also: DD_BASICPROC, DD_SNR

OK = false(size(Job.Nifti(SubjNr,:)));
if any(cellfun(@isempty,{Job.Nifti(SubjNr,:).Files}))
	return
end

% Launch a graphics window if we do a distributed job
dd_initcnode(Job)

% Get and save the (realigned) b0-mask
if strcmp(Job.DenoisingMenu.Str{Job.DenoisingMenu.Val}, 'none')
	Masks = dd_basicproc_getmask(Job, SubjNr, [], LogName);				% Create and print mask(s)
else
	Masks = dd_basicproc_getmask(Job, SubjNr);							% Create temporary mask(s): new denoised masks will be created later, do not print these
end

% Estimate the SNR and write some log info
for SeriesNr = 1:size(Job.Nifti,2)
	
	[b0Imgs DWImgs]			 = getb0imgs(Job, SubjNr, SeriesNr);	
	[SNR BGStd GhStd Spikes] = dd_snr(b0Imgs, Masks(SeriesNr), b0Imgs, 'Background noise (b0)');
	myspm_print(LogName, sprintf('S%d: Background noise (b0)', SeriesNr))
	FIDLog = fopen([LogName(1:end-2) 'tsv'], 'a');
	fprintf(FIDLog, 'S%d\tb0-SNR:\tSNR_SIG =\t%g\tSNR_SD =\t%g\n', SeriesNr, SNR);
	fprintf(FIDLog, 'S%d\tb0-GNR:\tGNR =\t%g\n', SeriesNr, mean(GhStd(:))/mean(BGStd(:)));
	fprintf(FIDLog, 'S%d\tb0-Spikes:\tn =\t%d\tN =\t%d\n', SeriesNr, Spikes);
	[SNR BGStd GhStd Spikes] = dd_snr(DWImgs, Masks(SeriesNr), DWImgs, 'Background noise (DWI)');
	myspm_print(LogName, sprintf('S%d: Background noise (DWI)', SeriesNr))
	fprintf(FIDLog, 'S%d\tDWI-SNR:\tSNR_SIG =\t%g\tSNR_SD =\t%g\n', SeriesNr, SNR);
	fprintf(FIDLog, 'S%d\tDWI-GNR:\tGNR =\t%g\n', SeriesNr, mean(GhStd(:))/mean(BGStd(:)));
	fprintf(FIDLog, 'S%d\tDWI-Spikes:\tn =\t%d\tN =\t%d\n', SeriesNr, Spikes);
	fclose(FIDLog);
	if ~strcmp(Job.DenoisingMenu.Str{Job.DenoisingMenu.Val}, 'none')
		delete(fullfile(Job.Nifti(SubjNr,SeriesNr).Path, 'mean*'))		% Clean-up the temporary mask(s)
	end
	
	% Return succesful
	OK(SeriesNr) = true;
	
end


function [b0Imgs DWImgs] = getb0imgs(Job, SubjNr, SeriesNr)
%
% Input:  dd_basicproc Job
% Output: filelist (SPM-style)

if isempty(Job.Nifti(SubjNr,SeriesNr).Files)
    warning('DIDI:OpenFile', ['Cannot find any images in: ' Job.Nifti(SubjNr,SeriesNr).Path])
	[b0Imgs DWImgs] = deal([]);
	return
end
for n = 1:numel(Job.Nifti(SubjNr,SeriesNr).Files)
	FList{n} = fullfile(Job.Nifti(SubjNr,SeriesNr).Path, Job.Nifti(SubjNr,SeriesNr).Files{n});
    D(n)	 = orderfields(dti_get_dtidata(FList{n}));
end
b0Sel  = [D.b]<=50;   % ==0;
b0Imgs = char(FList{b0Sel});
DWImgs = char(FList{~b0Sel});

if isempty(b0Imgs)
    warning('DIDI:OpenFile', ['Cannot find any b0-images in: ' fileparts(FList{1})])
end


function myspm_print(LogName, HdrTxt)

% Robust against closed figures
HG = spm_figure('FindWin', 'Graphics');
if ishandle(HG)
	HD = axes('Position', [0.05 0.95 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0.75, 0.5, datestr(now), 'Parent',HD)
	text(0, 0.5, HdrTxt, 'Parent',HD)
	spm_print(LogName)
end
