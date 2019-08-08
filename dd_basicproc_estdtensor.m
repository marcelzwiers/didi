function TDImgs = dd_basicproc_estdtensor(Job, SubjNr, SeriesNr, LogName)

% DD_BASICPROC_ESTDTENSOR is a key internal function of dd_basicproc that
% is made available externally to allow distributed computing.
%
% FUNCTION TDImgs = DD_BASICPROC_ESTDTENSOR(Job, SubjNr, SeriesNr, LogName)
%
% Marcel, 6-5-2011
%
% See also: DD_BASICPROC, DD_BASICPROC_REALIGNWARP

B0	   = 50;													% Threshold for selection of b0 images
TDImgs = '';
if strcmpi(Job.EstMenu.Str{Job.EstMenu.Val}, 'none') || isempty(Job.Nifti(SubjNr,SeriesNr).Files)
	return
end

HG		= dd_initcnode(Job);
rDWImgs = Job.Output(SubjNr,SeriesNr).TgtImgs;					% Fullpath cell array
Mask	= dd_basicproc_getmask(Job, SubjNr, SeriesNr);

% Avoid crashes of Volkmar's toolbox
if Job.ParallelBox.Val
	CmdLn = spm_get_defaults('cmdline');
	spm_get_defaults('cmdline', true)
end

% Limit the tensor estimation to data from the inner shell(s) (i.e. to higher SNR data)
if Job.MaxBValText.Str
	for n = 1:numel(rDWImgs)
		BVal(n) = getfield(dti_get_dtidata(rDWImgs{n}), 'b');
	end
	MaxBVal = str2double(Job.MaxBValText.Str);
	Sel     = BVal <= MaxBVal;
	rDWImgs = rDWImgs(Sel);
	fprintf('\n-> Limiting the tensor estimation to %i/%i DW images\n', sum(Sel), length(Sel));
	Job.Output(SubjNr,SeriesNr).WImgs = Job.Output(SubjNr,SeriesNr).WImgs(Sel(BVal > B0), :);
end

% Estimate the diffusion tensor
if strcmp(Job.EstMenu.Str{Job.EstMenu.Val}, 'GLM')
	if exist([Job.Nifti(SubjNr,SeriesNr).Path filesep 'SPM.mat'], 'file')
		delete([Job.Nifti(SubjNr,SeriesNr).Path filesep 'SPM.mat'])
	end
	ErrorVar = struct('erriid', 1);		% SPM5: = struct('errauto', struct('ltol',10, 'dtol',0, 'sep',0));
	EstMask	 = {''};					% SPM5: = Mask;	necessary for non-sphericity correction (see help dti_dt_regress)?
	PreFix	 = {'Dxx_' 'Dxy_' 'Dxz_' 'Dyy_' 'Dyz_' 'Dzz_'};
	dti_dt_regress(struct('srcimgs', {rDWImgs(:)}, 'errorvar', ErrorVar, ...
						  'dtorder',2, 'maskimg',{EstMask}, 'srcscaling','raw', ...
						  'swd',{{[Job.Nifti(SubjNr,SeriesNr).Path filesep]}}, 'spatsm',0))
elseif strcmp(Job.EstMenu.Str{Job.EstMenu.Val}, 'PATCH')
	PreFix = {'Dxx_PATCH_' 'Dxy_PATCH_' 'Dxz_PATCH_' 'Dyy_PATCH_' 'Dyz_PATCH_' 'Dzz_PATCH_'};
	DVols  = dd_wregress(rDWImgs, char(Job.Output(SubjNr,SeriesNr).WImgs), PreFix);		% DVols = xyzt volume
else
	error('Estimation option not recognized')
end
for n = 1:numel(PreFix)
	DImgs{n} = spm_file(rDWImgs{1}, 'prefix',PreFix{n});
end
if ~exist('DVols', 'var')
	DVols = spm_read_vols(spm_vol(char(DImgs)));
end

% Set the negative tensor eigenvalues to zero
NrVox  = dd_eigvalcorr(DImgs, Mask);
FIDLog = fopen([LogName(1:end-2) 'tsv'], 'a');
fprintf(FIDLog, 'S%g\tEigenvalue correction:\tNrCorr =\t%d\tNrVox =\t%d\n', SeriesNr, NrVox);
fclose(FIDLog);

% Display the tensor images
spm_check_registration(char(DImgs))
spm_orthviews('Interp', 0)					% No interpolation = better visibility of artefacts
spm_orthviews('context_menu','orientation',2)
spm_orthviews('XHairs','Off')
spm_orthviews('Window', [1 4 6], max(median(max(max(DVols(:,:,:,[1 4 6]))))) * [0 0.7])			% Use the same scaling for all diagonal tensor images
spm_orthviews('Window', [2 3 5], max(median(max(max(abs(DVols(:,:,:,[2 3 5])))))) * [-0.7 0.7])	% Use the same scaling for all off-diagonal tensor images
myspm_print(LogName, ['S' num2str(SeriesNr) ': Estimated tensor elements'])

% Compute the eigenvectors/values & radial diffusivity - TODO: plot eigenvectors & values
if Job.EigBox.Val
	dti_eig(struct('dtimg', {DImgs}, 'dteigopts', 'vl'))
	if Job.RDBox.Val
		% Load the eigenvalues, take the mean of the 2nd and 3rd value and save it as a eval23_*file
		Eval2  = strrep(DImgs{1},'Dxx_','eval2_');
		Eval3  = strrep(DImgs{1},'Dxx_','eval3_');
		Evals  = spm_vol(char(Eval2, Eval3));
		OutHdr = rmfield(Evals(1), 'pinfo');
		OutHdr.fname = strrep(Eval2,'eval2_','eval23_');
		spm_write_vol(OutHdr, mean(spm_read_vols(Evals),4));
	end
end

% Compute the tensor derivatives, export and display the results
TDVal = logical([Job.MDBox.Val Job.NormDBox.Val Job.NormABox.Val Job.FABox.Val Job.RABox.Val Job.ModeBox.Val Job.RDBox.Val Job.RDBox.Val]);	% If radial then also do axial.	NB: Match TDVal with DD_BASICPROC
if any(TDVal)
	TDOpts = 'dcnfvmar';
	TDImgs = char(strrep(DImgs{1},'Dxx_','ad_'),	strrep(DImgs{1},'Dxx_','nd_'), ...
				  strrep(DImgs{1},'Dxx_','na_'),	strrep(DImgs{1},'Dxx_','fa_'), ...
				  strrep(DImgs{1},'Dxx_','va_'),	strrep(DImgs{1},'Dxx_','mo_'), ...
				  strrep(DImgs{1},'Dxx_','eval1_'), strrep(DImgs{1},'Dxx_','eval23_'));
	TDOpts = TDOpts(TDVal);
	TDImgs = TDImgs(TDVal,:);
	if size(Job.Nifti,2)>1 && ~Job.FABox.Val && ~strcmp(Job.RealignMenu.Str{Job.RealignMenu.Val},'none')
		TDOpts = [TDOpts 'f'];		% Make sure we create an FA-image for hyper-alignment
	end
	dti_indices(struct('dtimg', {DImgs}, 'option', TDOpts))		% Catches TDOpts 'fvd'
	dd_indices2(struct('dtimg', {DImgs}, 'option', TDOpts))		% Catch TDOpts 'cnm'

	% Export the derivative images
	if Job.FSLBox.Val
		disp('-> Exporting masked tensor indices')
		for n = 1:size(TDImgs,1)
			[PName FName] = fileparts(TDImgs(n,:));
			system_dccn(sprintf(['source ~/.bashrc; fslchfiletype NIFTI_GZ %s %s; ' ...
				'fslmaths %s -mas %s/FDT_Data/nodif_brain_mask %s'], ...
				TDImgs(n,:), fullfile(PName,'FDT_Data',FName), ...
				fullfile(PName,'FDT_Data',FName), PName, fullfile(PName,'FDT_Data',FName)));
		end
	end

	% Display the derivate images
	spm_check_registration(TDImgs)
	spm_orthviews('context_menu','orientation',2)
	spm_orthviews('XHairs','Off')
	myspm_print(LogName, ['S' num2str(SeriesNr) ': Derived tensor indices'])

	% Display the derivate histograms
	NPlots	= size(TDImgs,1);				% We have maximal 8 derivatives, so plot up to 4 rows
	MaskVol = logical(spm_read_vols(spm_vol(char(Mask))));
	spm_figure('Clear','Graphics');
	for n = 1:NPlots
		ax	  = subplot(min(4,NPlots), ceil(NPlots/4), n, 'Parent',HG);
		TDVol = spm_read_vols(spm_vol(TDImgs(n,:)));
		hist(ax, TDVol(MaskVol), 100)
		title(spm_file(TDImgs(n,:),'filename'))
	end
	myspm_print(LogName, ['S' num2str(SeriesNr) ': Derived tensor indices'])
end

% Restore the original cmdline setting
if Job.ParallelBox.Val
	spm_get_defaults('cmdline', CmdLn)
end


%% ------------------- Auxillary Functions -------------------


function myspm_print(LogName, HdrTxt, FtrTxt)

% Robust against closed figures
HG = spm_figure('FindWin', 'Graphics');
if ishandle(HG)
	HD = axes('Position', [0.05 0.95 0.9 0.05], 'Visible','Off', 'Parent',HG);
	text(0.75, 0.5, datestr(now), 'Parent',HD)
	if nargin>1 && ~isempty(HdrTxt)
		text(0, 0.5, HdrTxt, 'Parent',HD)
	end
	if nargin>2 && ~isempty(FtrTxt)
		HD = axes('Position', [0.05 0 0.9 0.05], 'Visible','Off', 'Parent',HG);
		text(0, 0.5, FtrTxt, 'Parent',HD)
	end
	spm_print(LogName)
end