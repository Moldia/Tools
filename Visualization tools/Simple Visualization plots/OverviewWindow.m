% create an overview window
% useful for displaying plotted results
% Xiaoyan, 2015-6-30


hFig = gcf;
Figa = gca;
hLg = legend;
set(hLg,'Location','NorthEast','Color',[.6 .6 .6]);

hIm = Figa.Children;
hIm = hIm(end);



% create the scroll panel for interactive figure
hSP = imscrollpanel(hFig,hIm);
set(hSP,'Units','normalized',...
    'Position',[0 .05 1 .95]);

% create magnification box
hMagBox = immagbox(hFig,hIm);
pos = get(hMagBox,'Position');
set(hMagBox,'Position',[0 0 pos(3) pos(4)])
api = iptgetapi(hSP);
api.setVisibleLocation(1,1)

set(hFig,'toolbar','figure');

% create overview tool
imoverview(hIm)
api.setMagnification(api.findFitMag());