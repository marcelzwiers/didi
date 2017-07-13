function wTDImgs = dd_basicproc_normtdimgs(Job, SubjNr, LogName, TDImgs)

% DD_BASICPROC_NORMTDIMGS is a key internal function of dd_basicproc that
% is made available externally to allow distributed computing.
%
% FUNCTION TDImgs = DD_BASICPROC_NORMTDIMGS(Job, SubjNr, LogName)
%
% Marcel, 4-11-2014
%
% See also: DD_BASICPROC, DD_BASICPROC_ESTDTENSOR

% Defaults
wTDImgs = TDImgs;
TDVal	= logical([Job.MDBox.Val Job.NormDBox.Val Job.NormABox.Val Job.FABox.Val Job.RABox.Val Job.ModeBox.Val Job.RDBox.Val Job.RDBox.Val]);	% If radial then also do axial.	NB: Match TDVal with DD_BASICPROC
if strcmpi(Job.EstMenu.Str{Job.EstMenu.Val}, 'none') || ~any(TDVal) || any(cellfun(@isempty,{Job.Nifti(SubjNr,:).Files}))
	return
end
[HG HI]		  = dd_initcnode(Job);
NrSeries	  = size(Job.Nifti,2);
[Mask Meanb0] = dd_basicproc_getmask(Job, SubjNr);

% Clean-up the intermediate data files
Job = dd_basicproc_cleanup(Job, SubjNr, [], 'Denoise');

%--- Hyper-align the tensor derivative images
HyPar = cell(1, NrSeries);
if NrSeries>1 && ~strcmp(Job.HyperalignMenu.Str{Job.HyperalignMenu.Val},'none')
	
	disp('-> Hyper-aligning the tensor derivative images')

	% Create temporary masked FA volumes and an unbiased mean-image
	FAImgs = cell(NrSeries,1);
	for SeriesNr = 1:NrSeries
		[Pth FName Ext]	 = fileparts(deblank(TDImgs{SeriesNr}(1,:)));
		FAImgs{SeriesNr} = fullfile(Pth, ['fa_' FName(4:end) Ext]);
	end
	FAMask  = spm_vol(char(FAImgs));
	FAMask0 = FAMask;
	for n = 1:NrSeries
		FAVol			= spm_read_vols(FAMask(n));
		WVol			= spm_read_vols(spm_vol(Mask{n}));			% Todo: replace by MDVol? (bet-Mask is not very reliable)
		FAMask(n).fname = [tempname '.nii'];
		FAMask(n)		= spm_write_vol(FAMask(n), FAVol .* WVol);
	end
	spm_reslice(FAMask, struct('which',0, 'mean',1))
	MFAVol = spm_vol(spm_file(FAMask(1).fname, 'prefix','mean'));	% The images should already be in good register
	
	% Hyperalign the masked FA images
	switch Job.HyperalignMenu.Str{Job.HyperalignMenu.Val}
		case 'rigid-body'
			FAMask = spm_realign([MFAVol; FAMask], struct('quality',1, 'sep',2, 'rtm',1,'which',1));
		case 'affine'
			M = repmat(eye(4), [1 1 NrSeries]);
			for n = 1:NrSeries
				M(:,:,n)   = spm_affreg(MFAVol, FAMask(n), struct('regtype','rigid'));
				FAMask(n).mat = M(:,:,n) * FAMask(n).mat;			% Apply temporary realignment for reslicing to the mean
			end
			spm_reslice(FAMask, struct('which',0, 'mean',1))		% Save a new mean-FA volume for a two-pass procedure (as in spm_realign)
			MFAVol = spm_vol(spm_file(FAMask(1).fname, 'prefix','mean'));
			for n = 1:NrSeries
				FAMask(n).mat = M(:,:,n) \ FAMask(n).mat;			% Undo the temporary realignment
				M(:,:,n)	  = spm_affreg(MFAVol, FAMask(n), struct('regtype','rigid'), M(:,:,n));
				FAMask(n).mat = M(:,:,n) * FAMask(n).mat;
			end
		case 'non-linear'
			wTDImgs = prepend(wTDImgs, 'h');
			for SeriesNr = 1:NrSeries				% TODO: make this a two-pass procedure (as in affine)?
				HyPar{SeriesNr} = spm_file(FAMask0(SeriesNr).fname, 'suffix','_sn', 'ext','.mat');
				spm_normalise(MFAVol, FAMask(SeriesNr), HyPar{SeriesNr},[],[], struct('smoref',0, 'smosrc',0, 'regtype','rigid', 'reg',10));
				for n = 1:size(TDImgs{SeriesNr},1)
					spm_write_sn(TDImgs{SeriesNr}(n,:), HyPar{SeriesNr}, struct('bb',NaN(2,3), 'vox',NaN(1,3), 'prefix','h'))
				end
				spm_write_sn(Mask{SeriesNr}, HyPar{SeriesNr}, struct('bb',NaN(2,3), 'vox',NaN(1,3)))
			end
		otherwise
			error('Unknown hyper-align option: %s', Job.HyperalignMenu.Str{Job.HyperalignMenu.Val})
	end
	
	% Print and apply the linear / affine hyper-aligments to all derivative images
	disp('Hyperalignment:')
	for SeriesNr = 1:NrSeries
		if strcmp(Job.HyperalignMenu.Str{Job.HyperalignMenu.Val}, 'non-linear')
			HPar = load(HyPar{SeriesNr});
			HPar = spm_imatrix(HPar.Affine);
			fprintf('[%s]\n', sprintf('%.3f\t',HPar))
		else
			HPar = spm_imatrix(FAMask(SeriesNr).mat / FAMask0(SeriesNr).mat);
			fprintf('[%s]\n', sprintf('%.3f\t',HPar))
			[S LWarn] = mywarning('Off');						% Get rid of QFORM0 warning from affine transformations
			for n = 1:size(TDImgs{SeriesNr},1)
				spm_get_space(TDImgs{SeriesNr}(n,:), FAMask(SeriesNr).mat);
			end
			spm_get_space(Mask{SeriesNr}, FAMask(SeriesNr).mat);
			mywarning(S, LWarn)
		end
	end
	
	% Clean-up temporary FA & masked FA images (see DD_BASICPROC_ESTDTENSOR)
	delete(MFAVol.fname, FAMask.fname)
	if ~Job.FABox.Val
		delete(FAImgs{:})
	end
	
end

%--- Normalise the computed scalar results to the T1/T2 template and return the unwarped TDImgs
if Job.NormBox.Val
	
	disp('-> Normalizing the tensor derivative images to MNI space')
	
	% Estimate the normalization parameters
	spm('FigName', 'Normalizing DTI results', HI);
	T1Img = Job.Nifti(SubjNr).T1Text;					% Fullpath char array
	if Job.T1Box.Val && ~isempty(T1Img)					% We have a T1 image
		[Pth FName] = fileparts(T1Img);
		GMImg		= fullfile(Pth, ['wc1' FName '.' getfield(spm_get_defaults('images'),'format')]);	% NB: extension guessing not super robust
		if exist(GMImg, 'file')							% Be careful not to overwrite userfiles (i.e. temporarily move them elsewhere)
			TmpDir = tempname;
			mkdir(Pth, TmpDir)
			movefile(fullfile(Pth, ['wc1' FName '.*']), fullfile(Pth,TmpDir))
		end
		[Tag Job] = spm_jobman('harvest', spm_cfg_preproc8);
		Job.channel.vols = {T1Img};
		for n = 1:numel(Job.tissue)
			Job.tissue(n).native = [0 0];
			Job.tissue(n).warped = [0 0];
		end
		Job.tissue(1).warped = [1 0];					% Just save a GM-image to display, nothing else needed
		Job.warp.write		 = [1 1];					% Save the deformation fields
		Preproc{1}.spm.spatial.(Tag) = Job;
		SnPar = spm_jobman('run', Preproc);
	else												% Use the T2-template
		SnPar  = [Meanb0{1}(1:end-4) '_sn.mat'];
		NTempl = spm_vol(fullfile(spm('Dir'), 'toolbox', 'OldNorm', 'EPI.nii'));
		spm_normalise(NTempl, spm_vol(Meanb0{1}), SnPar);
		myspm_print(LogName, 'S*: Normalisation DWI -> EPI-template', Job.Nifti(SubjNr).Path)
	end
	
	% Normalize and display the TDImgs
	for SeriesNr = 1:NrSeries
		wTDImgs{SeriesNr} = my_write_sn2(TDImgs{SeriesNr}, SnPar, HyPar{SeriesNr});
		NTempl = fullfile(spm('Dir'), 'canonical', 'avg152T1.nii');
		spm_check_registration(char(NTempl, wTDImgs{SeriesNr}))
		if exist('GMImg', 'var')
			spm_orthviews('AddColouredImage', 1, GMImg, [1 0 0])		% [1 0 0] = 'red'
			spm_orthviews('Redraw')
		end
		myspm_print(LogName, ['S' num2str(SeriesNr) ': Derived tensor indices (warped to MNI-space)'])
	end
	
	% Clean up
	if exist('GMImg', 'var')
		delete([GMImg(1:end-3) '*'])
		if exist('TmpDir', 'var')						% Restore the original userfiles
			movefile(fullfile(Pth,TmpDir,'*'), Pth, 'f')
			rmdir(fullfile(Pth,TmpDir))
		end
	end
	
end


%% ------------------- Auxillary Functions -------------------


function wTDImgs = my_write_sn2(TDImgs, SnPar, HyPar)

% Combine the hyperaligmment with the normalization using SPM12's deformation utility
Def.comp{1}.sn2def.matname		= {HyPar};
Def.comp{1}.sn2def.vox			= NaN(1,3);
Def.comp{1}.sn2def.bb			= NaN(2,3);
if iscell(SnPar)
	Def.comp{2}.def				= SnPar{1}.fordef;	% ={fullfile(Pth, ['y_' FName '.nii'])}
else
	Def.comp{2}.sn2def.matname	= {SnPar};
	Def.comp{2}.sn2def.vox		= NaN(1,3);
	Def.comp{2}.sn2def.bb		= NaN(2,3);
end
if isempty(HyPar)
	Def.comp(1) = [];
end
Def.out{1}.pull.fnames			= cellstr(TDImgs);
Def.out{1}.pull.savedir.savesrc = 1;
Def.out{1}.pull.interp			= 4;
Def.out{1}.pull.mask			= 1;
Def.out{1}.pull.fwhm			= [0 0 0];
Job{1}.spm.util.defs			= Def;
Out		= spm_jobman('run', Job);
wTDImgs = char(Out{1}.warped);


function PO = prepend(PI, Pre)

PO = PI;
for n = 1:numel(PI)
	P = cell(size(PI{n},1),1);
	for m = 1:numel(P)
		[pth,nm,xt] = fileparts(char(PI{n}(m,:)));
		P{m}		= fullfile(pth, [Pre nm xt]);
	end
	PO{n} = char(P);
end


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