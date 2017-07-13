function [varargout] = df_reslice(ud,flags,msk,wip,MD2R)
%
% Reslices series of diffusion weighted images.
%
% FORMAT: [mask meanvol] = df_reslice(ud,flags,mask,wip,MD2R);
%
% INPUT:
% ud.P         - Handles of images to reslice.
% ud.beta      - Movement and distortion parameters
%                for images in ud.P.
% flags
% flags.mask   - Should we mask images so that only voxels for
%                which data exist in all scans are non-zero? 
%                Zero values occur when regions 'outside' the image 
%                are moved 'inside' the image during
% flags.hold   - Hold for interpolation (spm_sample_vol.m).
% flags.mean   - If set creates mean image.
% flags.which  - 0: Do not write any resliced images.
%                1: Write all resliced images.
%                2: Write mean image only.
%                3: Write all resliced images and mask image.
% wip          - Handle to function that reports on work in progress.
% MD2R		   - Transformation from diffusion to reference space (Marcel)
%_______________________________________________________________________
% Jesper Andersson 26/3-01

def_flags = struct('mask',         1,...
                   'hold',         1,...
                   'mean',         1,...
                   'which',        1);
if nargin < 1
   error('');
end

if nargin < 2 || isempty(flags)
   flags = def_flags;
else
   if ~isfield(flags,'mask'), flags.mask = def_flags.mask; end
   if ~isfield(flags,'hold'), flags.hold = def_flags.hold; end
   if ~isfield(flags,'mean'), flags.mean = def_flags.mean; end
   if ~isfield(flags,'which'), flags.which = def_flags.which; end
end

if nargin < 3 || isempty(msk);
   msk = 1:prod(ud.P(1).dim(1:3));
end

if nargin < 4
   wip = [];
end

if nargin < 5
   MD2R = eye(4);
end

[x,y,z] = ndgrid(1:ud.P(1).dim(1),1:ud.P(1).dim(2),1:ud.P(1).dim(3));
xyz = [x(:) y(:) z(:) ones(prod(ud.P(1).dim(1:3)),1)]';

%
% Evaluate largest common volume.
%
if flags.mask
   if ~isempty(wip), feval(wip,'DataMaskStart',length(ud.P)); end       
   fprintf('\n')
   for i=1:length(ud.P)
      b    = ud.beta((i-1)*(6+3*ud.P(i).dim(3))+1:i*(6+3*ud.P(i).dim(3)));
      txyz = df_transf_coord(ud.P(i),b,xyz(:,msk),[],MD2R); % Undistorted space + extra transformation
      msk  = msk(txyz(1,:)>=1 & txyz(1,:)<=ud.P(i).dim(1) &...
                 txyz(2,:)>=1 & txyz(2,:)<=ud.P(i).dim(2) &...
                 txyz(3,:)>=1 & txyz(3,:)<=ud.P(i).dim(3));
      if ~isempty(wip), feval(wip,'DataMaskUpdate',i); end
	  fprintf('.')
   end
   fprintf('\n')
   varargout{1} = msk;
end

%
% If we are not asked to write any resliced image
% our job is done here.
%
if ~flags.which && ~flags.mean
   return;
end

%
% Reslice images. Prefix resliced images with a 'u'.
%
if ~isempty(wip), feval(wip,'ResliceStart',length(ud.P)); end
if flags.mean 
   mvol = zeros(ud.P(1).dim(1:3)); 
   cvol = zeros(ud.P(1).dim(1:3)); 
end
for i=1:length(ud.P)
   b = ud.beta((i-1)*(6+3*ud.P(i).dim(3))+1:i*(6+3*ud.P(i).dim(3)));
   txyz = df_transf_coord(ud.P(i),b,xyz(:,msk),[],MD2R); % Undistorted space + extra transformation; Marcel
   f = spm_sample_vol(ud.P(i),txyz(1,:),txyz(2,:),txyz(3,:),flags.hold);
   ivol = zeros(ud.P(i).dim(1:3));
   ivol(msk) = f;
   if flags.mean
      mvol = mvol + ivol;
      cvol = cvol + (ivol > 0);
   end
   if flags.which == 1 || flags.which == 3
      p = ud.P(i);
      [path,fname,ext] = fileparts(p.fname);
      p.fname = fullfile(path,['u' fname ext]);
      p.descrip = [p.descrip ' Diffusion unwarped'];
	  % p.mat = MD2R * p.mat;		% Marcel
      p = spm_create_vol(p);
      spm_write_vol(p,ivol);
   end
   if ~isempty(wip), feval(wip,'ResliceUpdate',i,ud.P(i).fname); end
   fprintf('.')
end
if ~isempty(wip), feval(wip,'ResliceEnd'); end

%
% Write mean image if requested
%
if flags.mean
   mvol = mvol./cvol;
   varargout{2} = mvol;
   if flags.which ~= 0
      p = ud.P(1);
      [path,fname,ext] = fileparts(p.fname);
      p.fname = fullfile(path,['mean' fname ext]);
      p.descrip = [p.descrip ' Diffusion unwarped'];
	  % p.mat = MD2R * p.mat;		% Marcel
      p = spm_create_vol(p);
      spm_write_vol(p,mvol);
   end
end