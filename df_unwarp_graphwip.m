function df_unwarp_graphwip(mode,p1,p2,p3,p4)
%
% Will present information pertaining to the status of the
% iterative process of unwarping diffusion weighted images.
%
% FORMAT uw_graphwip(mode,p1,...)
%
% mode      - Verb specifying action.
% p1-p4     - Depends on mode.
%
%_______________________________________________________________________
% Jesper Andersson 26/3-01

persistent iter;
persistent resima;

if strcmp(mode(end-4:end), 'Start')
	HI = spm_figure('FindWin', 'Interactive');
	spm_figure('Clear', HI);
	spm('FigName', 'Affine realignment', HI);
end
switch mode
   case 'ItProcStart'
      iter = 0;
      resima = [];
   case 'IterStart'     % Report on start of new iteration.
      iter = p1;
   case 'SmoothStart'
      spm_progress_bar('Init',p1,'Creating temporary smooth images','Volumes smoothed');
   case 'SmoothUpdate'
      spm_progress_bar('Set',p1);
   case 'ImpMaskStart'    % Report on start of evaluation of importance of voxels
      spm_progress_bar('Init',p1,'Computing importance based mask','Slices finished');
   case 'ImpMaskUpdate'
      spm_progress_bar('Set',p1);
   case 'DataMaskStart'
      spm_progress_bar('Init',p1,...
        sprintf('Ascertaining data in all scans, iteration %d',iter-1),...
        'Scans checked');
   case 'DataMaskUpdate'
      spm_progress_bar('Set',p1);
   case 'DataMaskEnd'
      spm_progress_bar('Clear');
   case 'DerivativesStart'
      spm_progress_bar('Init',p1,...
        sprintf('Estimating partial derivatives, iteration %d',iter-1),'');
   case 'DerivativesUpdate'
      spm_progress_bar('Set',p1);
   case 'AtyStart'
      spm_progress_bar('Init',p1,...
        sprintf('Calculating A''y, iteration %d',iter-1),'Slices done');
   case 'AtyUpdate'
      spm_progress_bar('Set',p1);
   case 'IterEnd'  % End of one iteration
      if length(resima) < 3
         resima{length(resima)+1} = p4;
      else
         for i=1:2
            resima{i} = resima{i+1};
         end
         resima{3} = p4;
      end
      GraphDispIter(p1,p2,p3,resima);
   case 'ResliceStart'
      spm_progress_bar('Init',p1,'Reslicing images','Volumes done');
   case 'ResliceUpdate'
      spm_progress_bar('Set',p1);
   case 'ResliceEnd'
      spm_progress_bar('Clear');
   otherwise
      % Just ignore undefined mode commands.
end


%% END

function GraphDispIter(mss,beta,dim,resima)
%
% This is a utility function that Displays the result of 
% one iteration in the SPM graphics output window.
%
%_______________________________________________________________________
% Jesper Andersson 26/3-01


fg = spm_figure('FindWin','Graphics');
if isempty(fg)
	warning('Cannot find graphics window');
	return
end
spm_figure('Clear','Graphics');

%
% In top three panels we show profiles for
% (from left to right) shear, scale and
% translation. for first, middle and
% last diffusion gradient.
%

nscan = length(beta) / (6+3*dim(3));
beta = reshape(beta,6+3*dim(3),nscan);
ax = subplot(5,3,13, 'Parent', fg);
plot(ax, 1:dim(3),180/pi*asin(beta(7:6+dim(3),1)),...
		 1:dim(3),180/pi*asin(beta(7:6+dim(3),round(nscan/2))),...
		 1:dim(3),180/pi*asin(beta(7:6+dim(3),nscan)));
tmp = axis(ax); axis(ax, [1 dim(3) tmp(3:4)]);
title('Shear (degrees)');

ax = subplot(5,3,14, 'Parent', fg);
plot(ax, 1:dim(3),1+beta(6+dim(3)+1:6+2*dim(3),1),...
		 1:dim(3),1+beta(6+dim(3)+1:6+2*dim(3),round(nscan/2)),...
		 1:dim(3),1+beta(6+dim(3)+1:6+2*dim(3),nscan));
tmp = axis(ax); axis(ax, [1 dim(3) tmp(3:4)]);
title(ax, 'Scaling'); xlabel(ax, 'Slice #');

ax = subplot(5,3,15, 'Parent', fg);
plot(ax, 1:dim(3),beta(6+2*dim(3)+1:6+3*dim(3),1),...
		 1:dim(3),beta(6+2*dim(3)+1:6+3*dim(3),round(nscan/2)),...
		 1:dim(3),beta(6+2*dim(3)+1:6+3*dim(3),nscan));
tmp = axis(ax); axis(ax, [1 dim(3) tmp(3:4)]);
title(ax, 'Translation (mm)');

%
% Show images of objective function
% for history of iterations.
%

for slice = 1:3
	maxv = 0; minv = 1e6;
	for i = 1:length(resima)
		maxv = max(maxv,max(resima{i}.d(:,slice)));
		minv = min(maxv,min(resima{i}.d(:,slice)));
	end
	if minv==maxv
		minv = minv-eps;
	end
	for iter = 1:length(resima);
		ax = subplot(5,3,(3-slice)*3+iter, 'Parent', fg);
		imagesc(rot90(reshape(resima{iter}.d(:,slice),dim(1:2))), 'Parent',ax, [minv maxv]);
		axis(ax, 'image'); set(ax,'XTick',[],'YTick',[]);
		if iter == 1
			ylabel(ax, sprintf('Slice %d',resima{iter}.sl(slice)));
		end
		if slice == 3
			title(ax, sprintf('Iteration %d',length(mss)-length(resima)+iter-1));
		end
	end
end

%
% Show mss vs. iteration
%

ax = subplot(5,3,10, 'Parent', fg); plot(ax, 0:length(mss)-1,mss); title(ax, 'Objective function');
xlabel(ax, 'Iteration'); ylabel(ax, 'mss');

%
% Show rigid body parameters in lower right corner
%

ax = subplot(5,3,11, 'Parent', fg); plot(ax, beta(1:3,:)'); title(ax, 'Translation');
xlabel(ax, 'Scan #'); ylabel(ax, 'mm');
ax = subplot(5,3,12, 'Parent', fg); plot(ax, (180/pi)*beta(4:6,:)'); title(ax, 'Rotation');
xlabel(ax, 'Scan #'); ylabel(ax, 'degrees');
