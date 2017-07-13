function [txyz,jac] = df_transf_coord(P,beta,xyz,slmsk,MD2R,Inv)
%
% Returns transformed coordinate vector given rigid-body
% and slicewise eddy-current transforms in beta.
%
% Format:  txyz      = df_transf_coord(P,beta,[xyz],[slmsk])
% or      [txyz,jac] = df_transf_coord(P,beta,[xyz],[slmsk])
%                
%
% txyz       : nx3 matrix with transformed coordinates
% jac        : nx1 vector with voxel-wise Jacobian determinants.
%
% P          : Image volume struct (from spm_vol).
% beta       : 6+3*dim(3) parameter vector containing (in the 
%              following order) rigid-body movement parameters,
%              slicewise (xy) shears, slicewise (y) scaleings and
%              slicewise (y) translations.
% xyz        : 3xn matrix with original coordinates.
% slmsk      : dim(3)x1 vector specifying the number of voxels
%              for each slice that are contained in xyz.
% MD2R		 : Transformation from diffusion to undistorted 
%              reference (unweighted) space (Marcel)
% Inv		 : Computes the inverse transformation if true
%
%____________________________________________________________
% Jesper Andersson

if min(size(xyz))==0
   txyz = [];
   return
end

dim = P(1).dim(1:3);
if nargin < 3 || isempty(xyz)
   [x,y,z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
	xyz = [x(:) y(:) z(:) ones(prod(dim(1:3)),1)]';
end
% if nargin < 4 || isempty(slmsk)
%    if length(xyz) == prod(dim(1:3));
%       slmsk = cumsum(repmat(prod(dim(1:2)),1,dim(3)));
%    else
%       for i=1:dim(3)
%          slmsk(i) = length(find(xyz(3,:) == i));
%       end
%       slmsk = cumsum(slmsk);
%    end
% end
if nargin < 5 || isempty(MD2R)
   MD2R = eye(4);
end
if nargin<6 || isempty(Inv)
	Inv = false;
end
Inv = logical(Inv);

%
% Apply slice unique (eddy current) effects first.
%

%
% Edit by Marcel (correct, but slower):
%
if nargout > 1
   jac = zeros(1,size(xyz,2));
end
txyz	  = zeros(4,size(xyz,2));
T1		  = eye(4);									% Slicewise translation, shear & scaling
T2		  = spm_matrix(beta(1:6)');					% Rigid body transformation (of sampling grid)
iPM		  = eye(4) / P.mat;							% eye(4) = to ensure that the inverse is calculated outside the loop
T2iMD2RPM = T2 / MD2R * P.mat;						% MD2R = transformation of object/data
TT		  = P.mat \ T2 / MD2R * P.mat;

% 20110811 - Eelke
txyz = df_get_sampling_points(xyz, TT, iPM, T2iMD2RPM, beta', Inv);

return
%% The for-loop below is replaced by the df_get_sampling_points mex-file

for i = 1:dim(3)
	%
	% Get indexes to apply slice specific transform to.
	%
	if i==1
		indx = 1:slmsk(1);
	else
		indx = slmsk(i-1)+1:slmsk(i);
	end
	for n = indx
		ti = TT * xyz(:,n);							% Transformed sliceindex = P.mat \ T2 / MD2R * P.mat * [xyz(1:3,n); 1]
		ti = round(max([1 min([dim(3) ti(3)])]));	% 1 < ti < dim(3); NB: betas are linearly resampled in df_get_sampling_points instead (= better)
		T1(2,1)   = beta(6+ti);						% Slicewise shear
		T1(2,2)   = beta(6+ti+dim(3)) + 1;			% Slicewise scaling
		T1(2,4)   = beta(6+ti+2*dim(3));			% Slicewise translation
		if Inv
			txyz(:,n) = (iPM * T1 * T2iMD2RPM) \ xyz(:,n);	% = P.mat \ T1 * T2 / MD2R * P.mat
		else
			txyz(:,n) =  iPM * T1 * T2iMD2RPM  * xyz(:,n);	% = P.mat \ T1 * T2 / MD2R * P.mat
		end
	end
	
	% if nargout > 1								% This never happens
	% 	jac(indx) = 1 + beta(6+dim(3)+i);			% Correct?
	% end
end

txyz_ = df_get_sampling_points(xyz, TT, iPM, T2iMD2RPM, beta', Inv);
if sum(abs(txyz_(:)-txyz(:)))>1
	warning('df_get_sampling_points does not return the same sampling points')
end

return
%% Jesper's code (= wrong), suboptimally modified
%
% I agree with jesper the T1 should be applied before the T2, but believe
% that the order should be reversed when transforming *coordinates*!

if nargout > 1
   jac = zeros(1,size(xyz,2));
end
txyz = zeros(4,size(xyz,2));
T2 = spm_matrix(beta(1:6)');					% Rigid body transformation
for i=1:dim(3)
	ti = P.mat\T2/MD2R*P.mat * [0 0 i 1]';		% Transformed sliceindex
	ti = round(ti(3));
	if ti<1
		ti = 1;
	elseif ti>dim(3)
		ti = dim(3);
	end
	T1 = spm_matrix([0 beta(6+2*dim(3)+ti)]);	% Slicewise translation
	T1(2,1) = T1(2,1) + beta(6+ti);				% Slicewise shear
	T1(2,2) = T1(2,2) + beta(6+dim(3)+ti);		% Slicewise scaling
	%T = P.mat\T2*T1*P.mat;
	T = P.mat\T1*T2/MD2R*P.mat;
	%
	% Get indexes to apply slice specific transform to.
	%
	if i==1
		indx = 1:slmsk(1);
	else
		indx = slmsk(i-1)+1:slmsk(i);
	end
	if size(xyz,1) == 4
		txyz(:,indx) = T*xyz(:,indx);
	else
		txyz(:,indx) = T*[xyz(:,indx); ones(1,length(indx))];
	end
	if nargout > 1
		jac(indx) = 1 + beta(6+dim(3)+i);
	end
end
