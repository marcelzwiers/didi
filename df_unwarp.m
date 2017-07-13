function ud = df_unwarp(P,g,c,p,wip)
%
% Estimation of distortions of diffusion weighted images
%
% FORMAT: ud = df_unwarp([P],[g],[c],[p],[wip])
%
% INPUT:
% P             - Array of filenames or filehandles
% g             - length(P)x3 matrix with diffusion gradient vectors
% c             - Struct containing information on constraints.
% c.sptl        - Information on spatial constraints.
% c.sptl.apply  - Shall we use spatial constraints? (1/0)
% c.sptl.basset - Basis set ('cos'/'pol')
% c.sptl.order  - Order of spatial basis function
% c.grad        - Information on constraints across acquisitions.
% c.grad.apply  - Shall we use gradient constraints? (1/0)
% c.grad.shear  - Vector indicating which gradients we believe
%                 affect shearing. 
%                 c.grad.shear = [1] implies considering only 
%                 direct effects (self induction) of x-diffusion gradient.
%                 c.grad.shear = [1 2 3] implies that we
%                 model also cross terms, i.e. how a y- (2) and z-
%                 (3) diffusion gradient cause eddy currents in the
%                 "x-direction".
%                 c.grad.shear = [1 2 3 4 5 6 7 8 9] implies that a
%                 second order model is used, i.e. that squares and
%                 cross-products are used as well. The order of the
%                 terms is [x y z x^2 y^2 z^2 xy xz yz].
%                 c.grad.shear = [] implies that we do not constrain
%                 shears (though we may constrain scales and translations).
% c.grad.scale  - See c.grad.shear
% c.grad.trans  - See c.grad.shear
% p             - Struct containing "tuning" parameters.
% p.fwhm        - Width of smoothing kernel
% p.imask       - Information about "importance masking".
% p.imask.apply - Should "importance masking" be applied? (1/0)
% p.imask.type  -
% p.imask.prcnt -
% p.imask.vol	- additional mask (added by Marcel)
% p.difm        - Should we base R on diffusion tensor model? (1/0)
% p.orth        - Effects that should not be included in R.
%                 On might e.g. Orthogonalise DT-model design matrix 
%                 with respect to effects explicable by a first order 
%                 gradient constraint model. Depending on how one 
%                 selected ones directions, orthogonalisation may have 
%                 from very little (good directions) to very large
%                 (bad directions) effect.
% wip           - Handle to function that reports on work in progress.
%
% OUTPUT:
% ud            - Struct containing results from unwarping.
% ud.P          - Copy of P on input.
% ud.g          - Copy og g on input.
% ud.c          - Copy of c on input.
% ud.p          - copy of p on input.
% ud.beta       - Estimated distortion and movement parameters.
%                 On output beta is an n*(6+3*sl)_by_1 column
%                 vector where n is number of gradient directions and
%                 sl is number of slices. First comes the 6+3*sl parameters
%                 pertaining to the first scan, then the 6+3*sl pertaining
%                 to the second scan, etc. In each block of 6+3*sl
%                 parameters the 6 pertaining to movement comes first,
%                 followed by the sl parameters pertaining to shear, then
%                 scaling and finally translation.
% ud.mss        - Objective function for each iteration.
%_______________________________________________________________________
% Jesper Andersson 26/3-01

%
% Do a bit of input checking.
%
if nargin == 0
   P = spm_select(Inf,'image','Select diffusion weighted images');
end
if ~isfield(P(1),'dim'), P = spm_vol(P); end
[Dum1 Dum2 Ext] = fileparts(P(1).fname);
if nargin < 2, g = []; end
if size(g,2) == 3  % Prepare for second order model
   gg = [g g.^2 g(:,1).*g(:,2) g(:,1).*g(:,3) g(:,2).*g(:,3)];
end
if nargin < 3 
   c.sptl.apply = 0;
   c.grad.apply = 0;
end
if ~isfield(c,'sptl'), c.sptl.apply = 0; end
if ~isfield(c,'grad'), c.grad.apply = 0; end
if nargin < 4
   p.fwhm = 15;
   p.imask.apply = 1;
   p.imask.type = 2;
   p.imask.prcnt = 10;
   p.difm = 0;
   p.orth = [];
end
if nargin < 5 || isempty(wip)
   wip = @df_unwarp_textwip;
end
if isfield(p,'orth') && ~isempty(p.orth)
   if isfield(c,'grad') && isfield(c.grad,'apply') && c.grad.apply > 0
      error('Sorry, orthogonalisation and gradient constraints don''t work together');
   end
end

%
% Start building output struct
%

ud.P = P;
ud.g = g;
ud.c = c;
ud.p = p;

%
% Possibly work on smoothed images
%
if exist('p','var') && isfield(p,'fwhm') && ~isempty(p.fwhm) && p.fwhm > 0
   feval(wip,'SmoothStart',length(P));
   for i=1:length(P)
      feval(wip,'SmoothUpdate',i,P(i).fname);
      % [skrutt,tmpfname,Ext] = fileparts(tempname);
      % sfname(i,:) = fullfile(path,[tmpfname Ext]);
	  sfname(i,:) = [tempname Ext];
      spm_smooth(P(i).fname,sfname(i,:),p.fwhm);
   end
   P = spm_vol(sfname);
end
feval(wip,'SmoothEnd');


%
% Create matrices for spatial constraints
%

if c.sptl.apply == 1
	if strcmp(c.sptl.basset,'cos')
		B = spm_dctmtx(P(1).dim(3),c.sptl.order);
	elseif strcmp(c.sptl.basset,'pol')
		B = repmat([1:P(1).dim(3)]',1,c.sptl.order);
		for i=1:c.sptl.order
			B(:,i) = B(:,i).^(i-1);
			B(:,i) = B(:,i)./sum(B(:,i));
		end
		[B,skrutt] = qr(B,0);
		if B(1,1) < 0 B(:,1) = -1*B(:,1); end
	else
		error(sprintf('Invalid basis set %s',c.sptl.basset));
	end
	B = blkdiag(eye(6),B,B,B);
	BB = B;
	for i=2:length(P)
		BB = blkdiag(BB,B);
	end
end

%
% And matrix for "gradient" constraints
%

if c.grad.apply == 1
	if c.sptl.apply == 1
		sptl_order = c.sptl.order;
	else
		sptl_order = P(1).dim(3);
	end
	tmp = eye(6);     % Rigid part
	%
	% Shears
	%
	if ~isfield(c.grad,'shear') || isempty(c.grad.shear)
		tmp = blkdiag(tmp,eye(sptl_order));
	else
		tmp = [tmp; zeros(sptl_order,size(tmp,2))];
	end
	%
	% Scalings
	%
	if ~isfield(c.grad,'scale') || isempty(c.grad.scale)
		tmp = blkdiag(tmp,eye(sptl_order));
	else
		tmp = [tmp; zeros(sptl_order,size(tmp,2))];
	end
	%
	% Translations
	%
	if ~isfield(c.grad,'trans') || isempty(c.grad.trans)
		tmp = blkdiag(tmp,eye(sptl_order));
	else
		tmp = [tmp; zeros(sptl_order,size(tmp,2))];
	end
	
	affine = tmp;
	for i=2:length(P)
		affine = blkdiag(affine,tmp);
	end
	shear = [zeros(6,sptl_order); eye(sptl_order);...
		zeros(2*sptl_order,sptl_order)];
	scale = [zeros(6,sptl_order); zeros(sptl_order);...
		eye(sptl_order); zeros(sptl_order)];
	trans = [zeros(6,sptl_order); zeros(2*sptl_order,sptl_order);...
		eye(sptl_order)];
	G = affine;
	if isfield(c.grad,'shear') && ~isempty(c.grad.shear)
		G = [G kron(gg(:,c.grad.shear),shear)];
	end
	if isfield(c.grad,'scale') && ~isempty(c.grad.scale)
		G = [G kron(gg(:,c.grad.scale),scale)];
	end
	if isfield(c.grad,'trans') && ~isempty(c.grad.trans)
		G = [G kron(gg(:,c.grad.trans),trans)];
	end
end

%
% And residual forming matrix
%

if isfield(p,'difm') && p.difm == 1
	for i=1:length(P)
		X(i,:) = kron(g(i,:),g(i,:));
	end
	if isfield(p,'orth') && ~isempty(p.orth) && length(p.orth) == length(P)
		X = [(eye(length(P))-p.orth*pinv(p.orth)) * X ones(length(P),1)];
	end
	R = eye(length(P))-X*pinv(X);
else
	R = eye(length(P))-ones(length(P),1)*pinv(ones(length(P),1));
end

%
% Sampling grid
%

[x,y,z] = ndgrid(1:P(1).dim(1),1:P(1).dim(2),1:P(1).dim(3));
xyz = [x(:) y(:) z(:) ones(prod(P(1).dim(1:3)),1)]';
beta = zeros(length(P)*(6+3*P(1).dim(3)),1);


%
% Create mask based on "importance" of voxels.
%

if p.imask.apply ~= 0
	%
	% Always base mask on non-logged heavily smoothed data.
	%
	feval(wip,'ImpMaskStart',P(1).dim(3));
	mask_fname = [tempname Ext];
	mask_P = P(1);
	mask_P.fname = mask_fname;
	mask_P = spm_create_vol(mask_P);
	if isfield(p.imask, 'vol')	% Edit by Marcel: use extra brain-mask
		spm_write_vol(mask_P, spm_read_vols(P(1)) .* spm_dilate(p.imask.vol, ones([5,5,5])));
	else
		spm_write_vol(mask_P, spm_read_vols(P(1)));
	end
	spm_smooth(mask_fname, mask_fname, 20);
	mask_P = spm_vol(mask_fname);
	if p.imask.type == 1 % Clever but slow mask
		D = make_D(mask_P(1));
		importance = zeros(size(D,1),1);
		if c.sptl.apply == 1
			D = D*B;
		end
		D = D(:,[1 3:end]);
		DtD = D'*D;
		gc = det(inv(DtD));
		sz = prod(P(1).dim(1:2));
		for slice=1:P(1).dim(3)
			feval(wip,'ImpMaskUpdate',slice);
			for i=1:sz
				j = (slice-1)*sz+i;
				dtd = D(j,:)'*D(j,:);
				importance(j) = det(inv(DtD-dtd)) - gc;
			end
		end
		importance = importance ./ max(importance);
	elseif p.imask.type == 2 % Q&D mask
		[f,dx,dy,dz] = spm_sample_vol(mask_P(1),xyz(1,:),xyz(2,:),xyz(3,:),1);
		% D = make_D(mask_P(1));
		% importance = sqrt(dx.^2 + 2*dy.^2 + dz.^2);
		importance = sqrt(dx.^2 + dy.^2 + dz.^2);	% Edit by Marcel (y-factor 2 = typo?)
		clear f dx dy dz
	elseif p.imask.type == 3 % Really silly mask, for validation only.
		importance = spm_sample_vol(mask_P(1),xyz(1,:),xyz(2,:),xyz(3,:),0);
	end
	%
	% Keep most important percentage of voxels
	%
	[nhist,xhist] = hist(importance,1000);
	indx = find((cumsum(nhist)./sum(nhist)) > (1-p.imask.prcnt/100));
	threshold = xhist(indx(1));
	indx = find(importance > threshold);
	%
	% Connected component labelling discards FOV crap.
	%
	if p.imask.prcnt <= 20
		cc = spm_clusters(xyz(1:3,indx));
		[nhist,xhist] = hist(cc,max(cc));
		[skrutt,tcc] = max(nhist);
		imsk = indx(cc == tcc);
	else
		imsk = indx;
	end
end
feval(wip,'ImpMaskUpdate',P(1).dim(3));
feval(wip,'ImpMaskEnd');


%
% Here starts iterative search for warping parameters.
%

feval(wip,'ItProcStart');
for iter=1:10
	feval(wip,'IterStart',iter);
	%
	% Get mask where data exist for all scans.
	% Possibly combine with importance mask.
	%
	feval(wip,'DataMaskStart',length(P));
	if p.imask.apply == 0
		msk = 1:prod(P(1).dim(1:3));
	else
		msk = imsk;
	end
	for i=1:length(P)
		b = beta((i-1)*(6+3*P(i).dim(3))+1:i*(6+3*P(i).dim(3)));
		txyz = df_transf_coord(P(i),b,xyz(:,msk));
		indx = find(txyz(1,:)>=1 & txyz(1,:)<=P(i).dim(1) &...
					txyz(2,:)>=1 & txyz(2,:)<=P(i).dim(2) &...
					txyz(3,:)>=1 & txyz(3,:)<=P(i).dim(3));
		if isempty(indx)			% Runaway solution
			warning('Optimization failed (runaway solution): reverting to zero alignment')
			mss(iter) = mss(1);
			ud.mss	  = mss;
			ud.beta	  = zeros(length(P)*(6+3*P(1).dim(3)),1);
			return
		end
		msk = msk(indx);
		feval(wip,'DataMaskUpdate',i);
	end
	
	%
	% Let us use both a slicewise and a total mask.
	%
	for i=1:P(1).dim(3)
		slmsk(i) = length(find(xyz(3,msk) == i));	% = sum(xyz(3,msk)==i)
	end
	slmsk = cumsum(slmsk);
	
	%
	% Get partials of first scan.
	%
	feval(wip,'DerivativesStart',1);
	if exist('p','var') && isfield(p,'difm') && p.difm == 1
		[D,msk,slmsk] = make_D(P(1),beta(1:6+3*P(1).dim(3)),msk,slmsk,1);
	else
		D = make_D(P(1),beta(1:6+3*P(1).dim(3)),msk,slmsk,0);
	end
	feval(wip,'DerivativesUpdate',1);
	
	%
	% Get D'*D for single scan
	%
	
	DtD = D'*D;
	
	%
	% Introduce spatial constraints
	%
	
	if c.sptl.apply == 1
		DtD = B'*DtD*B;
		D = D*B;
	end
	
	%
	% And for entire problem
	%
	
	if c.grad.apply == 1
		AtA = G'*kron(R,DtD)*G;
	else
		AtA = kron(R,DtD);
	end
	
	%
	% Now calculate Aty
	%
	feval(wip,'AtyStart',P(1).dim(3));
	if ~exist('Temp', 'var')	% avoid updating
		Temp = find(diff([0 slmsk]));
	end
	% dispsl = [2 round(P(1).dim(3)/2) P(1).dim(3)-1];
	dispsl = [Temp(3) Temp(round(length(Temp)/2)) Temp(end-2)];		% Marcel
	dispevec.sl = dispsl;
	dispevec.d = zeros(prod(P(1).dim(1:2)),3);
	Aty = zeros(max(size(R).*size(DtD)),1);
	ss = 0;
	for sl=1:P(1).dim(3) % Do slicewise to save memory
		feval(wip,'AtyUpdate',sl);
		if sl==1
			indx = 1:slmsk(1);
		else
			indx = slmsk(sl-1)+1:slmsk(sl);
		end
		if isempty(indx)
			continue
		end
		evec = zeros(1,length(P)*length(indx));
		for i=1:length(P)
			if iter>1
				b = beta((i-1)*(6+3*P(i).dim(3))+1:i*(6+3*P(i).dim(3)));
				txyz = df_transf_coord(P(i),b,xyz(:,msk(indx)));
			else
				txyz = xyz(:,msk(indx));
			end
			f = spm_sample_vol(P(i),txyz(1,:),txyz(2,:),txyz(3,:),1);
			if exist('p','var') && isfield(p,'difm') && p.difm == 1
				f = log(f);
			end
			tmp = (f*D(indx,:))';
			for j=1:length(P)
				Aty((j-1)*size(D,2)+1:j*size(D,2)) =...
					Aty((j-1)*size(D,2)+1:j*size(D,2)) + R(j,i)*tmp;
				evec((j-1)*length(indx)+1:j*length(indx)) =...
					evec((j-1)*length(indx)+1:j*length(indx)) + R(j,i)*f;
			end
		end
		ss = ss + evec*evec';
		if any(sl == dispsl)
			tmp = sum(reshape(evec.*evec,length(indx),length(R)),2);
			tmpmsk = msk(indx) - (sl-1)*prod(P(1).dim(1:2));
			dispevec.d(tmpmsk, sl == dispsl) = tmp;
		end
	end
	if c.grad.apply == 1
		Aty = G'*Aty;
	end
	mss(iter) = ss / length(msk);
	
	%
	% Check for convergence. Dont bother with latest
	% estimates if reduction in mss too small.
	%

	if iter > 1 && (mss(iter-1)-mss(iter)) / mss(iter) < 1e-8
		break;
	end
	%
	% And get new beta
	%
	
	if c.grad.apply == 1
		beta0 = pinv(AtA)*Aty;
	else
		beta0 = kron(pinv(R),pinv(full(DtD)))*Aty;
	end
	
	%
	% Unpack beta if using constraints.
	%
	if c.grad.apply == 1
		beta0 = G*beta0;
	end
	if c.sptl.apply == 1
		beta0 = BB*beta0;
	end
	
	%
	% Update beta.
	%
	beta = beta - beta0;
	
	feval(wip,'IterEnd',mss,beta,P(1).dim(1:3),dispevec);
	fprintf(' iteration %g:%20s\n',iter-1,sprintf('mss = %f',mss(iter)));
end  % End of search for parameters

ud.mss = mss;
ud.beta = beta;


function [D,msk,slmsk] = make_D(P,beta,msk,slmsk,logflag)
%
% This is a utility function that builds a matrix of partial
% derivatives with respect to the parameters we are
% interested in.
%
% D is organised so that for every subject the six rigid transform
% parameters come first, followed by slicewise shears, slicewise
% scaleings and slicewise translations.
%_______________________________________________________________________
% Jesper Andersson 26/3-01

if nargin < 2 || isempty(beta)
   beta = zeros(6+3*P(1).dim(3));
end
if nargin < 3 || isempty(msk)
   msk = [1:prod(P(1).dim(1:3))];
end
if nargin < 4 || isempty(slmsk)
   slmsk = cumsum(repmat(prod(P(1).dim(1:2)),1,P(1).dim(3)));
end
if nargin < 5 || isempty(logflag)
   logflag = 0;
end

[x,y,z] = ndgrid(1:P(1).dim(1),1:P(1).dim(2),1:P(1).dim(3));
xyz = [x(:) y(:) z(:) ones(prod(P(1).dim(1:3)),1)]';

if all(beta == 0)
   txyz = xyz;
else
   txyz = df_transf_coord(P(1),beta(1:6+3*P(1).dim(3)),xyz);
end

tmpD = zeros(length(msk),8);
tiny = 0.00001;
if ~logflag
   [f,tmpD(:,1),tmpD(:,2),tmpD(:,3)] = spm_sample_vol(P(1),txyz(1,msk)',txyz(2,msk)',txyz(3,msk)',1);
else
   %
   % Note that we get the gradient in units of "per voxel".
   %
   for i=1:3
      ttxyz		 = txyz;
	  ttxyz(i,:) = txyz(i,:) + tiny*ones(1,length(txyz));
      tmp1		 = spm_sample_vol(P(1),ttxyz(1,msk),ttxyz(2,msk),ttxyz(3,msk),1);
      ttxyz(i,:) = txyz(i,:) - tiny*ones(1,length(txyz));
      tmp2		 = spm_sample_vol(P(1),ttxyz(1,msk),ttxyz(2,msk),ttxyz(3,msk),1);
      indx		 = find(tmp1 > 10 & tmp2 > 10);
      msk		 = msk(indx);
      tmpD		 = tmpD(indx,:);
      tmpD(:,i)  = (log(tmp1(indx)) - log(tmp2(indx)))' / (2*tiny);
   end
end

%
% In addition to the derivatives w.r.t. translations we want the
% derivatives with respect to rotation, shear and scaling.
% To that effect we define transformation matrices for these cases
%
tiny = 0.00001; 
for i=1:3
   tmp = zeros(1,3); tmp(i) = tiny;
   T{i} = spm_matrix([0 0 0 tmp]);
end
T{4} = eye(4); T{4}(2,1) = tiny;
T{5} = eye(4); T{5}(2,2) = 1+tiny;

%
% Rather than doing a numerical estimation of each derivative we
% simply rotate the gradient vectors. In this way
% we don't need separate code for logged and non-logged cases.
%
for i=1:5
   ttxyz = P(1).mat\T{i}*P(1).mat*txyz(:,msk);
   tmpD(:,3+i) = sum(((ttxyz(1:3,:)-txyz(1:3,msk))'.*tmpD(:,1:3)),2)/tiny;
end

%
% Get new slmsk in case logging has changed msk.
%
for i=1:P(1).dim(3)
   slmsk(i) = length(find(xyz(3,msk) == i));
end
slmsk = cumsum(slmsk);

% D = zeros(length(msk),6+3*P(1).dim(3));
D1 = zeros(length(msk),6);
%
% Rescale derivatives w.r.t. translation to "per mm"
%
for i=1:3
   tmp = zeros(1,3); tmp(i) = tiny;
   T = spm_matrix(tmp);
   ttxyz = P(1).mat\T*P(1).mat*txyz(:,msk);
   D1(:,i) =  sum(((ttxyz(1:3,:)-txyz(1:3,msk))'.*tmpD(:,1:3)),2)/tiny;
end

%
% Repackage rest of D
%
D1(:,4:6) = tmpD(:,4:6);

D = sparse(repmat([1:length(msk)]',6,1),reshape(repmat(1:6,length(msk),1),numel(D1),1),...
           reshape(D1,numel(D1),1),length(msk),6+3*P(1).dim(3));
clear D1;
        
col = zeros(length(msk),1);
for i=1:3
   for j=1:P(1).dim(3)
      if j==1
         indx = 1:slmsk(1);
      else
         indx = slmsk(j-1)+1:slmsk(j);
      end
      col(indx) = 6+(i-1)*P(1).dim(3)+j;
   end
   switch i
   case 1
      D = D + sparse(1:length(msk),col,tmpD(:,7),length(msk),6+3*P(1).dim(3));
   case 2
      D = D + sparse(1:length(msk),col,tmpD(:,8),length(msk),6+3*P(1).dim(3));
   case 3
      D = D + sparse(1:length(msk),col,D(:,2),length(msk),6+3*P(1).dim(3));
   end
end


