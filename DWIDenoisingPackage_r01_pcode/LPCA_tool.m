%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  LPCA demo tool
%
% Jose V. Manjon - jmanjon@fis.upv.es
% Pierrick Coupe - pierrick.coupe@labri.fr
%
% Copyright (C) 2008-2013 Jose V. Manjon and Pierrick Coupe 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all;
warning off;

addpath spm8

% parameters

rician=1;  % 1 for bias correction and 0 to disable it.
% number of logical CPUs (multithreading)
nbthread = maxNumCompThreads*2; % INTEL with old matlab
%nbthread=feature('numcores')*2; % INTEL with matlab 2013b
%nbthread=feature('numcores'); % AMD
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[filename, pathname] = uigetfile('*.nii;*.img', 'Select a nifty file/s','MultiSelect','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Proccess data ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~iscell(filename)) 
    temp=filename;
    filename=cell(1);
    filename{1}=temp;
end;
s=size(filename);

for i=1:s(2)
    
    disp(['Processing ',num2str(i),' of ',num2str(s(2))]);

    if(~isempty(pathname)) fname=[pathname,filename{i}];
    else fname=filename{i};
    end

    %read the data
    V=spm_vol(fname);
    ima=spm_read_vols(V);

    %filter the data        
    [fima] = DWIDenoisingLPCA(ima, rician, nbthread);

    % save result 
    ss=size(V);
    for ii=1:ss(1)
      V(ii).fname=[pathname,'lpca_',filename{i}];
      spm_write_vol(V(ii),fima(:,:,:,ii));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
