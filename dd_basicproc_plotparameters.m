function dd_basicproc_plotparameters(Params, Vols, LogName, HdrTxt, Headless)

% DD_BASICPROC_PLOTPARAMETERS is an internal function of dd_basicproc that is
% made available externally to allow distributed computing.
%
% Marcel, 26-8-2010
%
% See also: DD_BASICPROC, DD_BASICPROC_REALIGNPAR

% (Adapted from spm_realign.m)

% Use degrees instead of radials
Params(:,4:6) = Params(:,4:6)*180/pi;

disp('Transformation parameters:')
disp([{'x' 'y' 'z' 'pitch' 'rol' 'yaw'}; num2cell(round(1e4*Params(:,1:6))/1e4)])

for n = 1:numel(Vols)
	D(n) = dti_get_dtidata(Vols(n).fname);
end
b0s = find([D.b]<50);

fg = spm_figure('FindWin','Graphics');
if ishandle(fg)

    % display results
    % translation and rotation over time series
    %-------------------------------------------------------------------
    spm_figure('Clear',fg);
    ax = axes('Position',[0.1 0.7 0.8 0.2],'Parent',fg,'Visible','off');
	y = 1;
	text(0,y, 'DW Images','FontSize',16,'FontWeight','Bold','Parent',ax)
	y = y - 0.12;
	text(0,y, [fileparts(Vols(1).fname) ':'],'FontSize',10,'Interpreter','none','Parent',ax)
	y = y - 0.08;
	for i = 1:min([numel(Vols) 12])
		[Dum Name Ext] = fileparts(Vols(i).fname);
		text(0,y, sprintf('%-4.0f%s%s',i,Name,Ext),'FontSize',10,'Interpreter','none','Parent',ax)
		y = y - 0.08;
	end
	if numel(Vols) > 12
		text(0,y,'................ etc','FontSize',10,'Parent',ax)
	end
	
    ax = axes('Position',[0.1 0.4 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
    plot(Params(:,1:3),'Parent',ax)
	hold on
	plot(b0s(:),Params(b0s,1:3),'x','Parent',ax)
	xlim(ax, [1 size(Params,1)])
	set(ax, 'XTick', unique(round(get(ax,'XTick'))), 'YLimMode', 'Manual')
    set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold')
    set(get(ax,'Xlabel'),'String','image')
    set(get(ax,'Ylabel'),'String','mm')
	if ~Headless
		legend(ax, ['x';'y';'z'], 'Location','Best')	% Does not display well on headless machines
	end

    ax = axes('Position',[0.1 0.1 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
    plot(Params(:,4:6),'Parent',ax)
	hold on
	plot(b0s(:),Params(b0s,4:6),'x','Parent',ax)
	xlim(ax, [1 size(Params,1)])
	set(ax, 'XTick', unique(round(get(ax,'XTick'))), 'YLimMode', 'Manual')
    set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold')
    set(get(ax,'Xlabel'),'String','image')
    set(get(ax,'Ylabel'),'String','degrees')
	if ~Headless
		legend(ax, ['pitch';'roll ';'yaw  '], 'Location','Best')	% Does not display well on headless machines
	end
	
	HD = axes('Position', [0.05 0.95 0.9 0.05], 'Visible','Off', 'Parent',fg);
	text(0.75, 0.5, datestr(now), 'Parent',HD)
	if nargin>3 && ~isempty(HdrTxt)
		text(0, 0.5, HdrTxt, 'Parent',HD)
	end
	HD = axes('Position', [0.05 0 0.9 0.05], 'Visible','Off', 'Parent',fg);
	text(0, 0.5, LogName(1:end-3), 'Parent',HD)
	spm_print(LogName)
	
end
