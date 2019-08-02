function [HLeg,HLegBack,HAxBack] = resizeLegend(varargin)
% Goes to extreme lengths to create a shortened legend.
%
% @Synopsis
%   resizeLegend()
%     Does the above on the current axes for all plot elements.
%
%   resizeLegend(HAx)
%     Creates the legend on the given axes.
%
%   resizeLegend([HAx],'Name','Value')
%     Let's you adjust several other parameters.
%
%   [HLeg,HLegBack,HAxBack] = resizeLegend(...)
%     Returns the handles to the foreground legend and the background axes and
%     legend.
%
% @Input
%   HAx (optional)
%     A valid handle to an existing axes. If empty, the default is used.
%     DEFAULT: gca()
%
%   'LegendEntries',LegendEntries (Name/Value)
%     A cellstring with the display names, an array of handles for the plot
%     components which should be entered in the legend (as used by
%     LEGEND(entries)), simply 'show' (as used in LEGEND('show') or empty (in
%     which case the current legend is used as source for the entries).
%     DEFAULT: 'show' (all plot components)
%
%   'Scale',Scale (Name/Value)
%     A scalar factor >0 and <1 for which the legend entries should be shortened
%     (e.g. if a line is present in the legend (X1,X2) -> Scale*(X1,X2)). Values
%     below 0.3 will probably result in too short line segments.
%     DEFAULT: 0.5
%
%   'LegendProperties',LegendProperties (Name/Value)
%     A scalar struct with legend properties as passed to SET(HLegend,Props).
%     The struct must not contain the fields 'Box', 'Orientation', 'Parent',
%     'String' (use the input argument 'LegendEntries' instead) or 'Visible'. It
%     must also not contain a 'Position' and a 'Location' property
%     simultanously. If you specify the 'Position' property keep in mind that
%     the width and height of the legend will change and exact positioning
%     w.r.t. the righ/top corner will not be kept.
%     DEFAULT: struct()
%
% @Output
%   HLeg
%     A handle to the created legend.
%
%   HLegBack
%     A handle to the axes used as background for the legend. Its size is the
%     resized legend. Its 'Tag' property is set to 'LegendBackground'.
%
%   HAxBack
%     A handle to the axes used as background for the front-axes. Use this for
%     displaying a grid. Its 'Tag' property is set to 'AxesBackground'.
%
% @Remarks
%   - All input parameters can be passed as a struct (see
%     INPUTPARSER.STRUCTEXPAND).
%   - Minor & major grid have to be set on back-axis.
%   - For exact positioning, use HAxBack after calling RESIZELEGEND.
%
% @Dependencies
%   none;
%
% See also LEGEND
%
% ML-FEX ID:58914
%
% @Changelog
%   2016-08-29 (DJM): Beta release.
%   2016-09-08 (DJM):
%     - Now correclty treats legend entries of groups.
%     - Now correclty identifies the legend with mulitple plots in a figure.
%     - Does not change label font weights anymore.

ArgsOut = {'HLeg', 'HLegBack', 'HAxBack'};
nArgOut = nargout();

%% Check inputs
% Variable         : Type     : Acceptable Values       : Default
% -----------------:----------:-------------------------:--------
% HAx              : Optional : valid axes handle       : gca()
% LegendEntries    : Name/Val : nonempty                : 'show'
% Scale            : Name/Val : numeric scalar in (0,1) : 0.5
% LegendProperties : Name/Val : struct nonempty         : struct()

%Create INPUTPARSER with the passing of name/value arguments as struct enabled.
IP = inputParser();
IP.StructExpand = true;

%Validation functions.
validateHandleValid      = @(X)any([isempty(X) ~isempty(X)&ishandle(X)]);
validateStructScaNonempt = @(X)validateattributes(X,{'struct'},{'scalar','nonempty'});
validateNumScaInRange    = @(X)validateattributes(X,{'numeric'},{'scalar','>',0,'<',1});

%Define inputs.
IP.addOptional( 'HAx'             ,[] ,validateHandleValid);
IP.addParameter('LegendEntries'   ,'show');
IP.addParameter('Scale'           ,0.5     ,validateNumScaInRange);
IP.addParameter('LegendProperties',struct(),validateStructScaNonempt);

%Parse and validate all input arguments.
IP.parse(varargin{:});

%Assign optional & param/val arguments.
HAx        = IP.Results.HAx;
LegEntries = IP.Results.LegendEntries;
Scale      = IP.Results.Scale;
LegProps   = IP.Results.LegendProperties;

%Check axes & figure handle.
if isempty(HAx)
	HAx = gca();
end
HFig = HAx.Parent;

%Check if legend entries is cellstring.
if ~strcmpi(LegEntries,'show')
	if iscellstr(LegEntries)
		assert(all(cellfun(@(X) ischar(X),LegEntries)), ...
		'Legend entries must be a cell array of strings or an array of handles!');
	elseif isempty(LegEntries)
		HLeg = findobj(HFig.Children,'Type','legend');
		assert(ishghandle(HLeg), ['When legend entries are empty, the figure ' ...
			'must have a valid legend!']);
		LegEntries = HLeg.String;
		delete(HLeg);
	else
		assert(all(ishghandle(LegEntries)), ['Legend entries must be a cell ' ...
			'array of strings or an array of graphics object handles!']);
	end
end

%Check legend properties
for Field = {'Box','Orientation','Parent','String','Visible'}
	assert(~isfield(LegProps,Field{:}), ...
		'The ''%s'' property must not be specified in ''LegendProperties''!', ...
		Field{:});
end
assert(~(isfield(LegProps,'Position') & isfield(LegProps,'Location')), ...
	['The ''Position'' and the ''Location'' property must not both be ' ...
	'specified in ''LegendProperties''!']);
LegProps.Box     = 'off';
LegProps.Parent  = HFig;
LegProps.Visible = 'off';

%% Create background axes
HAxBack = copyobj(HAx,HFig);
uistack(HAxBack,'bottom');
cla(HAxBack);
set(HAxBack, 'Box','off', 'NextPlot','add', 'Tag','AxesBackground', ...
	'TickLength',[0 0], 'Title',text(), ...
	'XColor','none','YColor','none','ZColor','none', ...
	'XLabel',text(),'YLabel',text(),'ZLabel',text(), ...
	'XTickLabel','','YTickLabel','','ZTickLabel','');

%Make foreground axes transparent.
set(HAx,'Color','none','GridColor','none');
grid(HAx,'off');

%% Create legend
[HLeg,HIcons] = legend(LegEntries);
set(HLeg, LegProps);

%Get handles to the legend elements (use '-not' since ML creates a second, empty
%line segment for entries with no markers; thx to Sandro).
HLines   = findobj(HIcons,'Type','line', '-and','Marker'   ,'none', ...
                                         '-not','LineStyle','none');
HMarkers = findobj(HIcons,'Type','line', '-not','Marker'   ,'none', ...
                                         '-and','LineStyle','none');
HTexts   = findobj(HIcons,'Type','text');

%Check if HIcons contains a group of objects, in which case the tag of the
%lines/markers are empty.
if ~isempty(HLines)
    if isempty([HLines.Tag])
        set(HLines,{'Tag'},get(HTexts,'String'));
    end
end
if ~isempty(HMarkers)
    if isempty([HMarkers.Tag])
        set(HMarkers,{'Tag'},get(HTexts,'String'));
    end
end

%Store handles in the UserData.
HLeg.UserData.HAxes    = HAx;
HLeg.UserData.HLines   = HLines;
HLeg.UserData.HTexts   = HTexts;
HLeg.UserData.HMarkers = HMarkers;

%% Resize the legend elements
nEntries = length(HLeg.String);
MarkerTags = get(HMarkers,'Tag');
for i = 1:nEntries
	if i <= length(HLines)
		HLines(i).XData = Scale*HLines(i).XData;
		IsMarker = ismember(MarkerTags, HLines(i).Tag);
	else
		IsMarker = false(size(MarkerTags));
	end

	if i <= length(HMarkers)
		if any(IsMarker)
			HMarkers(i).XData = mean(HLines(IsMarker).XData);
		else
			HMarkers(i).XData = Scale*HMarkers(i).XData;
		end
	end

	if i <= length(HTexts)
		HTexts(i).Position(1) = Scale*HTexts(i).Position(1);
	end
end

%Compute the correct size of the resized legend.
Size = getLegSize(HLeg);

%% Create legend background
%Create legend background as third axes.
HLegBack = copyobj(HAxBack,HFig);
set(HLegBack, 'Box','off',...
              'GridColor','none',...
              'Units','normalized',...
              'Position',Size,...
              'Tag','LegendBackground',...
              'Visible','off',...
              'XColor','none',...
              'YColor','none',...
              'ZColor','none');
grid(HLegBack,'off');
uistack(HLegBack,'bottom');
uistack(HLegBack,'up',1);

%Store handle in the UserData.
HLeg.UserData.HLegBack = HLegBack;

%% Link properties
set([HAx;HAxBack;HLeg],'Units','normalized');

%Add callback functions for size recomputation if changed.
HFig.SizeChangedFcn = @sizeChangedCallback;
addlistener(HLeg,'Position','PostSet',@sizeChangedCallback);

linkprop([HAx;HAxBack],{'ActivePositionProperty',...
                        'GridLineStyle',...
                        'LineWidth',...
                        'MinorGridLineStyle',...
                        'OuterPosition',...
                        'Parent',...
                        'Position',...
                        'Units',...
                        'View',...
                        'XAxisLocation','YAxisLocation',...
                        'XDir','YDir','ZDir',...
                        'XLim','YLim','ZLim',...
                        'XMinorTick','YMinorTick','ZMinorTick',...
                        'XScale','YScale','ZScale',...
                        'XTick','YTick','ZTick'});
linkprop([HLeg;HLegBack], {'LineWidth','Units','Visible'});

%Make legend visible.
set([HLeg;HLegBack], 'Visible','on');

%% Check output
if nArgOut < length(ArgsOut)
	clearvars(ArgsOut{nArgOut+1:end});
end
end %RESIZELEGEND



%% SUBFUNCTION: getLegSize
function Size = getLegSize(HLeg)
% Computes the minimum size of the legend.
% @Input
%   HLeg (required): Handle to the legend.
%
% @Output
%   Size: New size in normalized units w.r.t. the figure as used in the
%         'Position' property of the legend.

%Set units to points to accomodate margin.
HFig = HLeg.Parent;
Units = get([HFig;HLeg],'Units');
set([HFig;HLeg],'Units','points');

ExtentsInPts = vertcat(HLeg.UserData.HTexts.Extent);
XSiz = inf;
if ~isempty(HLeg.UserData.HLines)
    XSiz = min([HLeg.UserData.HLines.XData, XSiz]);
end
if ~isempty(HLeg.UserData.HMarkers)
    XSiz = min([HLeg.UserData.HMarkers.XData, XSiz]);
end
if ~isfinite(XSiz)
    XSiz = 0;
end
Size = [XSiz ...                                %Min of lines & markers 
        ExtentsInPts(end,2) ...                 %Extent of last entry
        max(sum(ExtentsInPts(:,[1 3]),2)) ...   %Max extent of all entries
        sum(ExtentsInPts(1,[2 4]),2)];          %Max extent of first entry

%Convert reference frame to full axes & adjust for margin.
Size = Size.*HLeg.Position([3 4 3 4]);
Size(1:2) = HLeg.Position(1:2) + Size(1:2) ...
              - max([HLeg.UserData.HTexts.Margin])/2;

%Convert units to normalized w.r.t. the figure.
Size = Size./HLeg.Parent.Position([3 4 3 4]);
set([HFig;HLeg],{'Units'},Units);
end %GETLEGSIZE



%% SUBFUNCTION: sizeChangedCallback
function sizeChangedCallback(~,~)
% Recomputes the position of the legend and its background if they have changed.
	HFig = gcf;
	HLegs = findobj(HFig,'Type','legend');
	for i = 1:length(HLegs)
		if ~isempty(HLegs(i).UserData) && isequal(gca,HLegs(i).UserData.HAxes)
			Units = get(HLegs(i).UserData.HLegBack, 'Units');
			set(HLegs(i).UserData.HLegBack, 'Units','normalized');
			HLegs(i).UserData.HLegBack.Position = getLegSize(HLegs(i));
			set(HLegs(i).UserData.HLegBack, 'Units',Units);
		end
	end
end %SIZECHANGEDCALLBACK