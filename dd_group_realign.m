function QCPar = dd_group_realign(JobFile, RealignFile)

if nargin<1 || isempty(JobFile)
	JobFile = spm_select(1, 'mat');
end
if nargin<2 || isempty(RealignFile)
	RealignFile = spm_file(JobFile, 'suffix', '_realign', 'Ext','.tsv');
end

load(JobFile)
FID = fopen(RealignFile, 'w');
fprintf(FID, 'Subject\tSession\tFD\n');

for n = 1:size(Job.Nifti,1)
	for m = 1:size(Job.Nifti,2)
		
		QCFile = spm_file(JobFile, 'suffix',sprintf('_LogS%04d',n), 'Ext','.tsv');
		if ~exist(QCFile, 'file')
			continue
		end
		
		QCPar(n,m) = ParseQCPar(QCFile);
		FD		   = QCPar(n,m).Trans_X + QCPar(n,m).Trans_Y + QCPar(n,m).Trans_Z + QCPar(n,m).Rot_X + QCPar(n,m).Rot_Y + QCPar(n,m).Rot_Z;
		fprintf(FID, '%s\t%d\t%f\n', Job.Nifti(n,m).Path, m, FD);
		
	end
end

fclose(FID);


function QCPar = ParseQCPar(QCFile)

% Read the QA-parameters from the *.tsv log-files
FID  = fopen(QCFile);
if FID==-1
	QCPar = struct();
	warning('Failed to open: %s', QCFile)
	return
end

Line = fgetl(FID);
while ischar(Line)
	
	ID = textscan(Line, '%s', 'Delimiter','\t');
	switch ID{1}{2}
		case 'b0-mask:'							% dd_basicproc_getmask.m
			QCPar.TBVol			= cell2mat(textscan(Line, '%*s%*s%*s%f', 'Delimiter','\t'));
		case 'b0-SNR:'							% dd_snr.m
			% QCPar.SNR			= cell2mat(textscan(Line, '%*s%*s%*s%f%*s%*f', 'Delimiter','\t'));		% SNR_SIG works better with the LPCA filter. TODO: use LPCA noise estimate?
			QCPar.SNR			= cell2mat(textscan(Line, '%*s%*s%*s%*f%*s%f', 'Delimiter','\t'));		% SNR_SD works better for SONATA data
		case 'b0-GNR:'							% dd_snr.m
			QCPar.GNR			= cell2mat(textscan(Line, '%*s%*s%*s%f', 'Delimiter','\t'));
		case 'b0-Spikes:'						% dd_snr.m
			QCPar.Spikes		= 0;	% First image sometimes has elevated BG noise. cell2mat(textscan(Line, '%*s%*s%*s%f%*s%*f', 'Delimiter','\t'));
		case 'DWI-SNR:'							% dd_snr.m
			QCPar.SNR_DWI		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*s%f', 'Delimiter','\t'));
		case 'DWI-GNR:'							% dd_snr.m
			QCPar.GNR_DWI		= cell2mat(textscan(Line, '%*s%*s%*s%f', 'Delimiter','\t'));
		case 'DWI-Spikes:'						% dd_snr.m
			QCPar.Spikes		= (cell2mat(textscan(Line, '%*s%*s%*s%f%*s%*f', 'Delimiter','\t')) + QCPar.Spikes) * 1000 / cell2mat(textscan(Line, '%*s%*s%*s%*f%*s%f', 'Delimiter','\t'));
		case 'mean(DWI)-mask:'
			QCPar.Volume		= cell2mat(textscan(Line, '%*s%*s%*s%f', 'Delimiter','\t'));
		case 'Realignment (raw DWI):'			% dd_basicproc_realignpar.m
			QCPar.Trans_X_raw	= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%f%*f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans_Y_raw	= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans_Z_raw	= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot_X_raw		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%*f%f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot_Y_raw		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%*f%*f%f%*f', 'Delimiter','\t'));
			QCPar.Rot_Z_raw		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%*f%*f%*f%f', 'Delimiter','\t'));
		case 'PATCH:'							% dd_patch.m
			QCPar.Outl_Vox		= cell2mat(textscan(Line, '%*s%*s%*s%f%*s%*f%*s%*f', 'Delimiter','\t'));
			QCPar.Outl_Slc		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*s%f%*s%*f', 'Delimiter','\t'));
		case 'Coregistration DWI->T1:'			% dd_basicproc_realign.m
			QCPar.Trans2T1_X	= cell2mat(textscan(Line, '%*s%*s%*s%f%*f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans2T1_Y	= cell2mat(textscan(Line, '%*s%*s%*s%*f%f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans2T1_Z	= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot2T1_X		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot2T1_Y		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%f%*f', 'Delimiter','\t'));
			QCPar.Rot2T1_Z		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%f', 'Delimiter','\t'));
		case 'Coregistration DWI->b0:'			% dd_basicproc_realign.m
			QCPar.Trans2b0_X	= cell2mat(textscan(Line, '%*s%*s%*s%f%*f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans2b0_Y	= cell2mat(textscan(Line, '%*s%*s%*s%*f%f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans2b0_Z	= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot2b0_X		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot2b0_Y		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%f%*f', 'Delimiter','\t'));
			QCPar.Rot2b0_Z		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%f', 'Delimiter','\t'));
		case 'Realignment (patched):'			% dd_basicproc_realign.m
			QCPar.Trans_X		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%f%*f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans_Y		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%f%*f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Trans_Z		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%f%*f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot_X			= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%*f%f%*f%*f', 'Delimiter','\t'));
			QCPar.Rot_Y			= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%*f%*f%f%*f', 'Delimiter','\t'));
			QCPar.Rot_Z			= cell2mat(textscan(Line, '%*s%*s%*s%*f%*f%*f%*f%*f%*f%*s%*f%*f%*f%*f%*f%f', 'Delimiter','\t'));
		case 'abs(PEUnwarpDef):'				% dd_basicproc_warp.m
			QCPar.PEDefMean		= cell2mat(textscan(Line, '%*s%*s%*s%f%*s%*f', 'Delimiter','\t'));
			QCPar.PEDefMax		= cell2mat(textscan(Line, '%*s%*s%*s%*f%*s%f', 'Delimiter','\t'));
		case 'Eigenvalue correction:'			% dd_basicproc_estdtensor.m
			QCPar.NrEigCorr		= cell2mat(textscan(Line, '%*s%*s%*s%f%*s%*f', 'Delimiter','\t'));
		otherwise
			warning('Unknown option: %s in %s', ID{1}{2}, QCFile)
	end
	Line = fgetl(FID);
	
end
fclose(FID);
