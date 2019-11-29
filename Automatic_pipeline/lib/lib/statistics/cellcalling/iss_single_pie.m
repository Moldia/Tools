function iss_single_pie(iCell, o)
% iss_single_pie(iCell, o)

% find classes to collapse
CollapseMe = zeros(length(o.ClassNames),1);
Colors = zeros(length(o.ClassNames),3);
DisplayName = o.ClassNames;
for i=1:size(o.ClassCollapse,1)
    ClassList = o.ClassCollapse{i,1};
    for j=1:length(ClassList)
        MyClasses = strmatch(ClassList{j}, o.ClassNames);
        if length(MyClasses)==0; continue; end
        CollapseMe(MyClasses)=i;
        Colors(MyClasses,:) = repmat(o.ClassCollapse{i,3},length(MyClasses),1);
        DisplayName(MyClasses) = o.ClassCollapse(i,2);
    end
end

nColorWheel = sum(CollapseMe==0);

Colors0 = hsv(ceil(nColorWheel*1.2));
Colors(~CollapseMe,:) = Colors0(1:nColorWheel,:); % last is zero

figure; 
set(gcf, 'Color', 'w');
set(gca, 'color', 'w');
hold on

pMy = o.pCellClass(iCell,:);

% merge three subtypes of PC.CA1
pMy(1) = sum(pMy(1:3));
pMy(2:3) = 0;

WorthShowing = find(pMy>o.MinPieProb);
if ~isempty(WorthShowing)
    clf;
    h = pie(pMy(WorthShowing), DisplayName(WorthShowing));
    coordCell = o.CellYX(iCell,:);
    
    for i=1:length(h)/2
        hpatch = (i*2-1);
        %         Index = find(strcmp(h(i*2).String, NickNames), 1);
        set(h(hpatch), 'FaceColor', Colors(WorthShowing(i),:));
        % size based on number of reads
        set(h(hpatch), 'Xdata', get(h(hpatch), 'Xdata')*o.PieSize*sum(o.pSpotCell(:,iCell)) + coordCell(2));
        set(h(hpatch), 'Ydata', get(h(hpatch), 'Ydata')*o.PieSize*sum(o.pSpotCell(:,iCell)) + coordCell(1));
        
        set(h(hpatch), 'EdgeAlpha', 1, 'LineWidth', .1, 'EdgeColor', [.5 .5 .5]);
        
        htext = i*2;
        set(h(htext), 'Position', h(htext).Position*o.PieSize*sum(o.pSpotCell(:,iCell)) + [fliplr(coordCell) 0],...
            'FontSize', 12);
    end
    set(gca, 'XLim', coordCell(2)+[-50 50], 'YLim', coordCell(1)+[-50 50]);
    
end

axis equal
axis off
title(['Cell' num2str(iCell)]);
drawnow;

end
