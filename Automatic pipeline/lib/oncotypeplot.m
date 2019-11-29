function oncotypeplot(Risk, Polygon, windowsize, scale, img)
RG = {'low', 'g'; 'intermediate', 'y'; 'high', 'r'};

%% prepare UI
S.fh = figure('units','pixels','position',[200 200 900 800],...
    'menubar','figure',...
    'name','OncotypeDX',...
    'numbertitle','off',...
    'resize','on');
S.x = windowsize;  % For plotting.
S.ax = axes('unit','normalized',...
    'position',[.1 .1 .8 .8]);

% initial plot
S.f = cell(0,1);
if length(img)>1
    axes(S.ax);
    imshow(img)
else
    set(S.ax,'YDir','reverse');
end

plotoverlay(length(windowsize))

S.sl = uicontrol('style','slide',...
    'unit','normalized',...
    'position',[.1 .05 .8 .02],...
    'min',windowsize(1),'max',windowsize(end),'val',windowsize(end),...
    'sliderstep',[(windowsize(2)-windowsize(1))/range(windowsize) 1],...
    'callback',{@sl_call,S});


    function sl_call(varargin)
        % delete the previously plotted patches first
        cellfun(@(v) delete(v),S.f)

        % Callback for the slider.
        [h,S] = varargin{[1,3]};  % calling handle and data structure.
        v = windowsize==round(get(h,'value')/(windowsize(2)-windowsize(1)))*(windowsize(2)-windowsize(1));

        plotoverlay(v)
    end

    function plotoverlay(v)
        hold on;
        for j = 1:3
            temp = strcmp(Risk{v},RG{j,1});
            temp = Polygon{v}(:,:,temp);
            for k = 1:size(temp,3)
                temph = patch(temp(:,1,k),temp(:,2,k),RG{j,2},...
                    'facealpha',0.3,'edgealpha','0.1');
                S.f = [S.f, {temph}];
            end
        end
        title(['window size ' num2str(windowsize(v)/scale)])
        
    end
end
