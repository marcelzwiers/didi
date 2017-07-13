function [NiList T1Text] = dd_basicproc_convert2nifti(Job, SubjNr)

% DD_BASICPROC_CONVERT2NIFTI is an internal function of dd_basicproc that is
% made available externally to allow distributed computing.
%
% Marcel, 26-8-2010
%
% See also: DD_BASICPROC

% DICOM => Nifti file conversion

% Sort and convert DICOM files
CWD	= pwd;
cd(Job.DICOM(SubjNr).Path)
for n = 1:numel(Job.DICOM(SubjNr).Files)
	DcmList{n} = fullfile(Job.DICOM(SubjNr).Path, Job.DICOM(SubjNr).Files{n});
end
DcmList = char(DcmList);
DcmHdrs = spm_dicom_headers(DcmList);
Order = [];
for n = 1:numel(DcmHdrs)
	Order(n,:) = [DcmHdrs{n}.SeriesNumber   DcmHdrs{n}.AcquisitionNumber ...
				  DcmHdrs{n}.InstanceNumber DcmHdrs{n}.EchoNumbers];
end
[Dum Idx] = sortrows(Order);
DcmHdrs	  = DcmHdrs(Idx);
DcmList	  = DcmList(Idx,:);
try
	Out	= spm_dicom_convert(DcmHdrs,'all','flat',getfield(spm_get_defaults('images'),'format'));	% Works with spm8 and higher
catch
	Out = spm_dicom_convert_mz2(DcmHdrs);				% Robust against anonymous files
end
NiList = char(Out.files);
M1 = spm_get_space(NiList(1,:));
if any(abs(diag(M1(1:3,1:3))./eig(M1(1:3,1:3))) < 0.5)	% Reorient to xyz ordening if non-diagonal (e.g. sagittal scans)
	fprintf('--> Reorienting converted nifti files in: %s\n', fileparts(NiList(1,:)))
	switch spm_file(NiList(1,:),'ext')
		case 'nii'
			OutputType = 'NIFTI';
		otherwise
			OutputType = 'NIFTI_PAIR';
	end
	for NiFile = NiList'
		system(['source ~/.bashrc; export FSLOUTPUTTYPE=' OutputType '; fslreorient2std ' NiFile' ' ' NiFile']);
	end
end
if Job.T1Box.Val && ~isempty(Job.DICOM(SubjNr).T1Text)
	T1Ext = spm_file(Job.DICOM(SubjNr).T1Text(1,:),'ext');
	if strcmp(T1Ext, 'nii') || strcmp(T1Ext, 'img')
		T1Text = Job.DICOM(SubjNr).T1Text;
	else
		cd(fileparts(Job.DICOM(SubjNr).T1Text(1,:)))
		T1Hdrs = spm_dicom_headers(Job.DICOM(SubjNr).T1Text);
		if strcmp(spm('ver'), 'SPM5')
			T1Text = spm_dicom_convert_mz(T1Hdrs);
		else
			try
				Out = spm_dicom_convert(T1Hdrs,'all','flat',getfield(spm_get_defaults('images'),'format'));	% Works with spm8 and higher
			catch
				Out = spm_dicom_convert_mz2(T1Hdrs);
			end
            T1Text = char(Out.files);
		end
	end
else
	T1Text = [];
end

% Clean up old files and attach DW info to the Nifti files
NrImPerVol  = numel(DcmHdrs) / size(NiList,1);		% =1 for mosaic-dicom
[BVal BVec] = dicominfodti(DcmList(1:NrImPerVol:end,:));
if isempty(BVal)
	BVal = load(fullfile(Job.DICOM(SubjNr).Path, 'bval.txt'));
	BVec = load(fullfile(Job.DICOM(SubjNr).Path, 'bvec.txt'));	% Assume these are in world coordinates
	Sz	 = size(Bvec);
	if ~any(Sz==3) || numel(Sz~=2)
		error('The b-vector file has wrong size [%d x %d] instead of [n x 3] or [3 x n]', Sz)
	elseif Sz(1)==3
		BVec = BVec';
		BVal = BVal';
	end
end
[D DimFlip] = dd_rotategradients;
for n = 1:size(NiList,1)
	dti_get_dtidata(NiList(n,:), struct('mat',spm_get_space(NiList(n,:)), ...
		'b',BVal(n), 'g',DimFlip.*BVec(n,:)));		% Include Volkmar's x-flip error! (see dd_rotategradients)
	[Pth FName] = fileparts(NiList(n,:));
	delete(fullfile(Pth, ['*_' FName '*']), fullfile(Pth, ['*' FName '_*']))
end

cd(CWD)
