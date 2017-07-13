function gendiffvecfile(Sets)

% FUNCTION gendifvecfile(<Sets>)
%
% Generates a DiffusionVectors.txt in the current directory using the
% evenly distributed and optimally directions used in CAMINO. See:
%
% INPUT
%   Sets - Array that specifies the nr of DW directions
%          (default = [12 18 30 60 61 64 256])
%
% Cook PA, Symms M, Boulby PA and Alexander DC, Optimal acquisition orders
% of diffusion-weighted MRI measurements, Journal of Magnetic Resonance
% Imaging, 25(5), 1051-1058, 2007.
%
% Multiple unweighted (b0) images are interleaved in this sequence.
%
% Marcel, 27-2-2008.

if nargin<1 || isempty(Sets)
    Sets = [12 18 30 60 61 64 256];
end

FID = fopen('DiffusionVectors.txt', 'wt');
fprintf(FID, ['# ----------------------------------------------------------------------\n' ...
              '#\n' ...
              '# File: %%CustomerSeq%%\\DiffusionVectors.txt\n' ...
              '# Version: 1.0\n' ...
              '# Author: Marcel Zwiers\n' ...
              '# Date: ' datestr(now) '\n' ...
              '#\n' ...
              '# Lang: Input file for SBBDiffusion\n' ...
              '#\n' ...
              '# Descrip: This file contains external vector files\n' ...
              '# that are evenly distributed and optimally ordered. See:\n' ...
              '#\n' ...
              '# Cook PA, Symms M, Boulby PA and Alexander DC, Optimal\n' ...
              '# acquisition orders of diffusion-weighted MRI measurements,\n' ...
              '# Journal of Magnetic Resonance Imaging, 25(5), 1051-1058, 2007.\n' ...
              '#\n' ...
              '# Multiple unweighted (b0) images are interleaved in this sequence\n' ...
              '#\n' ...
              '# ----------------------------------------------------------------------\n']);

% 12 DW-directions
if any(Sets==12)
    Dirs = load('/home/mrphys/marzwi/Documents/DiffusionDirections/Elec012_ordered.txt');
    printvecs(FID, Dirs, 13)
end

% 18 DW-directions
if any(Sets==18)
    Dirs = load('/home/mrphys/marzwi/Documents/DiffusionDirections/Elec018_ordered.txt');
    printvecs(FID, Dirs, 10)
end

% 30 DW-directions
if any(Sets==30)
    Dirs = load('/home/mrphys/marzwi/Documents/DiffusionDirections/Elec030_ordered.txt');
    printvecs(FID, Dirs, 10)
end

% 60 DW-directions
if any(Sets==60)
    Dirs = load('/home/mrphys/marzwi/Documents/DiffusionDirections/Elec060_ordered.txt');
    printvecs(FID, Dirs, 10)
end

% 61 DW-directions
if any(Sets==61)
    Dirs = load('/home/mrphys/marzwi/Documents/DiffusionDirections/Elec061_ordered.txt');
    printvecs(FID, Dirs, 10)
end

% 64 DW-directions + 1 b0 (for Elena?)
if any(Sets==64)
    Dirs = load('/home/mrphys/marzwi/Documents/DiffusionDirections/PointSets/Elec064.txt');
    printvecs(FID, Dirs, 100)
end

% 256 DW-directions
if any(Sets==256)
    Dirs = load('/home/mrphys/marzwi/Documents/DiffusionDirections/Elec256_ordered.txt');
	Dirs = Dirs';
	Dirs = [Dirs(1); Dirs(4:end)'];
    printvecs(FID, Dirs, 12)
end

fclose(FID);
system('unix2dos DiffusionVectors.txt');


function printvecs(FID, Dirs, Rem)

% Every Rem there will be a b0 vector

Dirs(1) = [];  % Contains the number of directions
Dirs    = reshape(Dirs, [3 numel(Dirs)/3]);

NVecs = size(Dirs,2);
Nb0   = ceil(NVecs/(Rem-1));
fprintf(FID, ['\n# %g DW directions + %g b-zeros (interleaved)\n' ...
              '[directions=%g]\n' ...
              'CoordinateSystem = xyz\n' ...
              'Normalisation = none\n'], NVecs, Nb0, NVecs+Nb0);
B0Cntr = 0;
for n = 0:NVecs+Nb0-1
    if rem(n,Rem)==0
        fprintf(FID, 'Vector[%g] = ( 0, 0, 0 )\n', n);
        B0Cntr = B0Cntr + 1;
    else
        BVec = Dirs(:,n+1-B0Cntr);
        fprintf(FID, 'Vector[%g] = ( %g, %g, %g )\n', n, BVec/norm(BVec));
    end
end
