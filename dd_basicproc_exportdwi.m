function OK = dd_basicproc_exportdwi(Job, SubjNrs, SeriesNrs)

% DD_BASICPROC_EXPORTDWI is a key internal function of dd_basicproc that
% is made available externally to allow distributed computing.
%
% FUNCTION DD_BASICPROC_EXPORTDWI(Job, SubjNr, SeriesNr)
%
% Marcel, 1-3-2011
%
% See also: DD_BASICPROC, DD_BASICPROC_REALIGNWARP

if ~Job.CaminoBox.Val && ~Job.FSLBox.Val
	OK = true;
	return
end
if nargin<3 || isempty(SeriesNrs)
	SeriesNrs = 1:size(Job.Nifti,2);
end
if nargin<2 || isempty(SubjNrs)
	SubjNrs = 1:size(Job.Nifti,1);
end

for SubjNr = SubjNrs
	for SeriesNr = SeriesNrs

		% Get the list of DW images that need to be exported
		DWImgs = Job.Output(SubjNr,SeriesNr).TgtImgs;
		if isempty(DWImgs)
			continue
		end

		% Recompute the DW info in voxel coordinates (FSL and Camino do not consider the volume orientation)
		M1	  = spm_get_space(DWImgs{1});					% M1 = voxel-to-world mapping
		VoxSz = spm_imatrix(M1);
		VoxSz = VoxSz(7:9);
		for n = 1:numel(DWImgs)
			D(n) = orderfields(dti_get_dtidata(DWImgs{n}));
		end
		[D DimFlip] = dd_rotategradients(D, inv(M1));		% The dimflip brings Volkmar's vector in world-space (see dd_rotategradients)
		BVecs = diag(sign(VoxSz).*DimFlip) * vertcat(D.g)';	% Volkmar-space => image-space => voxel-space
		BVals = [D.b];

		[PName FName] = fileparts(DWImgs{1});
		% Export to Camino
		if Job.CaminoBox.Val
			CamDir = [PName filesep 'Camino_Data' filesep];
			if ~exist(CamDir, 'dir')
				mkdir(CamDir);
			else
				delete(fullfile(CamDir, '*'))
			end
			Vol = permute(spm_read_vols(spm_vol(char(DWImgs))), [4 1 2 3]);
			FID = fopen([CamDir FName '.Bfloat'], 'w', 'b');
			fwrite(FID, Vol, 'float');
			fclose(FID);
			disp('-> Exporting bvals and bvecs')
			FID = fopen([CamDir 'scheme'], 'w');
			fprintf(FID, '# Scheme file created by dd_basicproc\nVERSION: BVECTOR\n');
			for n = 1:length(BVals)
				fprintf(FID, '%f\t%f\t%f\t%G\n', BVecs(:,n)', BVals(n)*1e6);
			end
			fclose(FID);
		end

		% Export to FSL (NB: the fsl-calls can easily be replaced by spm-calls if needed)
		% NB: For export of masked FA-image see DD_BASICPROC_ESTDTENSOR (FA not computed yet)
		if Job.FSLBox.Val
			Mask   = dd_basicproc_getmask(Job, SubjNr, SeriesNr);
			FSLDir = fullfile(PName,'FDT_Data',filesep);
			if ~exist(FSLDir, 'dir')
				mkdir(FSLDir);
			else
				delete(fullfile(FSLDir, '*'))
			end
			disp(['-> Merging DWI files to FSL''s 4D-nifti format in ' FSLDir])
			system(['source ~/.bashrc; fslmerge -t ' FSLDir 'data' sprintf(' "%s"', DWImgs{:})]);
			disp('-> Exporting Mask')
			system(['source ~/.bashrc; fslchfiletype NIFTI_GZ ' char(Mask) ' ' FSLDir 'nodif_brain_mask']);
			disp('-> Exporting Masked (skull-stripped) mean b0-image')
			system(['source ~/.bashrc; fslchfiletype NIFTI_GZ ' strrep(char(Mask), '_mask', '') ' ' FSLDir 'nodif_brain']);
			disp('-> Exporting bvals and bvecs')
			save([FSLDir 'bvecs'], 'BVecs', '-ASCII');
			save([FSLDir 'bvals'], 'BVals', '-ASCII');
		end

	end
end

OK = true;
