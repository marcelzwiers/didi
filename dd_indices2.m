function varargout = dd_indices2(bch)

% Compute various invariants of diffusion tensor
%
% Batch processing
% FORMAT res = dd_indices2(bch)
% ======
%/*\begin{itemize}*/
%/*\item*/ Input argument:
%   bch struct containing fields
%/*\begin{description}*/
%/*\item[*/       .dtimg/*]*/  cell array of filenames of DT images
%/*\item[*/       .data/*]*/   Xdim-by-Ydim-by-Zdim-by-6 array of tensors
%               Only one of .files and .data need to be given.
%/*\item[*/       .option/*]*/ character array with a combination of
%/*\begin{description}*/
%/*\item[*/               c/*]*/ Tensor Norm
%/*\item[*/               n/*]*/ Anisotropy Norm
%/*\item[*/               m/*]*/ Anisotropy Mode
%/*\end{description}\end{description}*/
%/*\item*/ Output argument:
%   res struct containing fields (if requested)
%/*\begin{description}*/
%/*\item[*/       .nd/*]*/ ND output
%/*\item[*/       .na/*]*/ NA output
%/*\item[*/       .mo/*]*/ MO output
%/*\end{description}*/
%   where output is a filename if input is in files and a
%   Xdim-by-Ydim-by-Zdim array if input is a data array.
%/*\end{itemize}*/

rev = '$Revision$';
funcname = 'Anisotropy/Diffusivity';

% function preliminaries
Finter = spm_figure('FindWin', 'Interactive');
spm_input('!DeleteInputObj');
spm('FigName', funcname, Finter);
spm('FnBanner', mfilename, rev);

% function code starts here

if ~isfield(bch,'data')
	IDT = spm_vol(char(bch.dtimg));
	X = spm_read_vols(IDT);
	[p n e] = fileparts(IDT(1).fname);
	if ~isempty(regexp(n, '^((dt[1-6])|(D[x-z][x-z]))_.*'))
		n = n(5:end);
	end
else
	X = bch.data;
end
AD = mean(X(:,:,:,[1 4 6]), 4);
A = X;								% The anisotropic (deviatoric) part of the tensor
A(:,:,:,1) = A(:,:,:,1) - AD;
A(:,:,:,4) = A(:,:,:,4) - AD;
A(:,:,:,6) = A(:,:,:,6) - AD;
if any(bch.option=='n') || any(bch.option=='m')
	NA = sqrt((sum(A(:,:,:,[1 4 6]).^2,4) + 2*sum(A(:,:,:,[2 3 5]).^2,4)));
end

if any(bch.option=='c')
	ND = sqrt((sum(X(:,:,:,[1 4 6]).^2,4) + 2*sum(X(:,:,:,[2 3 5]).^2,4)));
	if ~isfield(bch, 'data')
		res.nd = {fullfile(p, ['nd_' n e])};
		IND = IDT(1);
		IND.fname = res.nd{1};
		IND.descrip = 'Tensor Norm';
		IND.pinfo = [1; 0; 0];
		IND.dt = [spm_type('float32') spm_platform('bigend')];
		IND = spm_create_vol(IND);
		spm_write_vol(IND, ND);
	else
		res.nd = ND;
	end
end
if any(bch.option=='n')
	if ~isfield(bch, 'data')
		res.na = {fullfile(p, ['na_' n e])};
		INA = IDT(1);
		INA.fname = res.na{1};
		INA.descrip = 'Anisotropy Norm';
		INA.pinfo = [1; 0; 0];
		INA.dt = [spm_type('float32') spm_platform('bigend')];
		INA = spm_create_vol(INA);
		spm_write_vol(INA,NA);
	else
		res.na = NA;
	end
end
if any(bch.option=='m')
	MO = 3 * sqrt(6) * (A(:,:,:,1) .* (A(:,:,:,4).*A(:,:,:,6) - A(:,:,:,5).*A(:,:,:,5)) - ...
						A(:,:,:,2) .* (A(:,:,:,2).*A(:,:,:,6) - A(:,:,:,5).*A(:,:,:,3)) + ...
						A(:,:,:,3) .* (A(:,:,:,2).*A(:,:,:,5) - A(:,:,:,4).*A(:,:,:,3))) ./ (NA.^3);
	if ~isfield(bch, 'data')
		res.mo = {fullfile(p, ['mo_' n e])};
		IMO = IDT(1);
		IMO.fname = res.mo{1};
		IMO.descrip = 'Anisotropy Mode';
		IMO.pinfo = [1; 0; 0];
		IMO.dt = [spm_type('float32') spm_platform('bigend')];
		IMO = spm_create_vol(IMO);
		spm_write_vol(IMO,MO);
	else
		res.mo = MO;
	end
end

spm_input('!DeleteInputObj');
if nargout>0
	varargout{1} = res;
else
	if isfield(bch,'data')
		fprintf('extracted information saved into workspace variable ''res''\n');
		assignin('base', 'res', res);
		disp(res);
	end
end
