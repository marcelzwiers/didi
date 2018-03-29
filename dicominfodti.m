function [BVal BVec BMtx Mosaic] = dicominfodti(DCMFiles, ASCIIFiles, CoordSys)

% [BVal BVec BMtx Mosaic] = dicominfodti(DCMFiles, ASCIIFiles, CoordSys)
%
% Reads b-values and gradient directions from Siemens dicom headers
%
% ASCII files containing BVal (names bval) and BVec (named bvec) will be written
% in the same directory when ASCIIFiles is true (default=false) or in the
% location specified by ASCIIFiles (character array).
%
% CoordSys determines the coordinate system in which the gradient vectors
% are represented. CoordSys can be RAH/RAS (default) or LAH/LAS (experimental).
%
% The results are displayed (including the erroneous b-vectors) when no
% output arguments are given.
%
% NB: Siemens private CSA fields are used. The BVecs are calculated from
% the b-matrix and *NOT* read directly from the header as this field is
% unreliable (at least in some (older) VB software versions).
%
% Marcel, 2-12-2008.

% From: http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:DICOM_for_DWI_and_DTI
%
% The recommended tags to use in DICOM are as follows:
% 0018 9075 CS 1 Diffusion Directionality
% 0018 9076 SQ 1 Diffusion Gradient Direction Sequence
% 0018 9087 FD 1 Diffusion b-value
% 0018 9089 FD 3 Diffusion Gradient Orientation
% 0018 9117 SQ 1 MR Diffusion Sequence
% 0018 9147 CS 1 Diffusion Anisotropy Type
%
% Private vendor: Siemens
% 0019;000C;SIEMENS MR HEADER  ;B_value                         ;1;IS;1
% 0019;000D;SIEMENS MR HEADER  ;DiffusionDirectionality         ;1;CS;1
% 0019;000E;SIEMENS MR HEADER  ;DiffusionGradientDirection      ;1;FD;3
% 0019;000F;SIEMENS MR HEADER  ;GradientMode                    ;1;SH;1
% 0019;0027;SIEMENS MR HEADER  ;B_matrix                        ;1;FD;6

if nargin<1 || isempty(DCMFiles)
	DCMFiles = spm_select(Inf, 'any', 'Select DICOM files', {''}, pwd, '(?i).ima$|(?i).dcm$');
end
if nargin<2 || isempty(ASCIIFiles)
    ASCIIFiles = true;
end
if nargin<3 || isempty(CoordSys)
	% See: http://eeg.sourceforge.net/mri_orientation_notes.html and http://trackvis.org/blog/forum/diffusion-toolkit-usage/siemens-gradient-orientation-and-dtk/
	CoordSys = 'RAS';
end

DCMHdrs = spm_dicom_headers(DCMFiles);
for n = 1:length(DCMHdrs)
    if ~checkfields(DCMHdrs{n},'ImageType','CSAImageHeaderInfo') ||...
            isempty(read_acquisitionmatrixtext(DCMHdrs{n})) ||...
            isempty(read_nrofimagesinmosaic(DCMHdrs{n}))
        Mosaic(n) = false;
        warning('DIDI:DicomInfo', 'Extracting DTI information was tested for mosaic images only...\n%s', DCMHdrs{n}.Filename)
    else
        Mosaic(n) = true;
    end
    try
        BVal(n) = str2num(char(cellstr(get_numaris4_val(DCMHdrs{n}.CSAImageHeaderInfo, 'B_value'))));
    catch
        warning('This file does not seem to contain DTI-info\n%s', DCMHdrs{n}.Filename)
		if nargout, [BVal BVec BMtx Mosaic] = deal([]); end
		return
    end
    if BVal(n)==0
        BVec(n,:) = [0 0 0];
    else
        try
            BMtx(n,:) = str2num(char(cellstr(get_numaris4_val(DCMHdrs{n}.CSAImageHeaderInfo, 'B_matrix')))); % BMtx = [bxx bxy bxz byy byz bzz] (?)
            [Y I]     = max(BMtx(n, [1 4 6]));      % Avoid using a zero cross-term (credits to Phil Cook)
            switch I
                case 1
                    BSign = sign(BMtx(n,1:3));      % Assume bx is positive
                case 2
                    BSign = sign(BMtx(n,[2 4 5]));  % Assume by is positive
                case 3
                    BSign = sign(BMtx(n,[3 5 6]));  % Assume bz is positive
            end
            BVec(n,:) = BSign .* sqrt(BMtx(n, [1 4 6])/BVal(n));
            % BVec(n,:) = BMtx(n,1:3)/sqrt(BVal(n)*BMtx(n,1));      % This is how Siemens does it
        catch
            BMtx      = [];
            BVec(n,:) = str2num(char(cellstr((get_numaris4_val(DCMHdrs{n}.CSAImageHeaderInfo, 'DiffusionGradientDirection'))))); % This is not reliable!
			if n==1
				disp('WARNING: The DICOM b-matrix field is not present. The b-vector information may be unreliable...')
			end
        end
    end
end
if strcmpi(CoordSys, 'RAH') || strcmpi(CoordSys, 'RAS')
    % Go to the RAS-coordinate system (DICOM = LPS)
    BVec(:, [1 2]) = -BVec(:, [1 2]);
elseif strcmpi(CoordSys, 'LAH') || strcmpi(CoordSys, 'LAS')
    % Go to the LAS-coordinate system (DICOM = LPS)
    BVec(:, 1) = -BVec(:, 1);
end
BVal = BVal(:);

if ASCIIFiles
	if ischar(ASCIIFiles)
		PathStr = ASCIIFiles;
	else
		PathStr = fileparts(DCMFiles(1,:));
	end
    save(fullfile(PathStr, 'bval.txt'), 'BVal', '-ASCII')
    save(fullfile(PathStr, 'bvec.txt'), 'BVec', '-ASCII')
end

if nargout==0    % Display the b-values and vectors from the header
    if isempty(BMtx)
        fprintf('Bval\tBvec\t\t\tNorm\n')
    else
        fprintf('Bval\tBVec(Mtx)\t\tNorm\t\tBvec\t\t\tNorm\n')
    end
    for n = 1:length(BVal)
        if isempty(BMtx)
            fprintf('%d\t%.3f\t%.3f\t%.3f\t%.3f\n', BVal(n), BVec(n,:), norm(BVec(n,:)))
        else
            if BVal(n)==0
                BVec2 = [0 0 0];
            else
				try
					BVec2 = str2num([DCMHdrs{n}.CSAImageHeaderInfo(22).item(1:3).val]); % This is not reliable!
					BVec2(1:2) = -BVec2(1:2);  % Go to our RAS-coordinate system
				catch
					BVec2 = NaN(1,3);
				end
            end
            fprintf('%d\t%.3f\t%.3f\t%.3f\t%.3f\t\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                 BVal(n), BVec(n,:), norm(BVec(n,:)), BVec2, norm(BVec2))
        end
    end
    clear all
end

%% ----------------------------------------------

function Ok = checkfields(DCMHdr, varargin)
Ok = true;
for n = 1:(nargin-1)
    if ~isfield(DCMHdr, varargin{n})
        Ok = false;
        break
    end
end

%% ----------------------------------------------

function dim = read_acquisitionmatrixtext(DCMHdr)
str = DCMHdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'AcquisitionMatrixText');
dim = sscanf(val','%d*%d')';
if length(dim)==1,
	dim = sscanf(val','%dp*%d')';
end;
if isempty(dim), dim=[]; end


%% ----------------------------------------------

function n = read_nrofimagesinmosaic(DCMHdr)
str = DCMHdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'NumberOfImagesInMosaic');
n   = sscanf(val','%d');
if isempty(n), n=[]; end;


%% ----------------------------------------------

function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str),
	if strcmp(deblank(str(i).name),name),
		for j=1:str(i).nitems,
			if  str(i).item(j).xx(1),
				val = {val{:} str(i).item(j).val};
			end;
		end;
		break;
	end;
end;
val = char(val{:});
