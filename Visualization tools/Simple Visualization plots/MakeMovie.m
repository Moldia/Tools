
% open the image .fig file first
hgload('Plotted_FullRes.fig');

% prepare
writerObj = VideoWriter('mBrain_short_low.mp4','MPEG-4');   % H.264 encoding
% writerObj = VideoWriter('mBrain_short.avi');
% set(writerObj,'Quality',100);
open(writerObj);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

zoom out;
FrameRegion = [127 35 1124 800]; 

% initial frames
for i = 1:20
    frame = getframe(gca,FrameRegion);
    writeVideo(writerObj,frame);
end

% zoom in
for k = 1:60
    zoom(1.02)
    for i = 1:2
        drawnow;
        frame = getframe(gca,FrameRegion);
        writeVideo(writerObj,frame);
    end
end

% stand for a little while
for i = 1:10
    frame = getframe(gca,FrameRegion);
    writeVideo(writerObj,frame);
end

% move to certain region
target = [6470 12800];  % starting x-,y- point of ROI
viewregion_x = get(gca,'XLim');
viewregion_y = get(gca,'YLim');

% left
while viewregion_x(1)>target(1)
    viewregion_x = viewregion_x - range(viewregion_x)*0.005;
    axis([viewregion_x,viewregion_y]);
    
    for i = 1:1
        drawnow;
        frame = getframe(gca,FrameRegion);
        
        writeVideo(writerObj,frame);
    end
end

% stand
for i = 1:5
    frame = getframe(gca,FrameRegion);
    writeVideo(writerObj,frame);
end

% up
while viewregion_y(1)<target(2)
    viewregion_y = viewregion_y + range(viewregion_y)*0.005;  % the image was flipped in Y direction to get a better representation (i.e. get(gca,'YDir')='normal' instead of 'reverse' as in most cases related to images)
    axis([viewregion_x,viewregion_y]);
    
    for i = 1:1
        drawnow;
        frame = getframe(gca,FrameRegion);
        writeVideo(writerObj,frame);
    end
end

% stand
for i = 1:20
    frame = getframe(gca,FrameRegion);    
    writeVideo(writerObj,frame);
end

% zoom in 
for k = 1:20
    zoom(1.1)
    for i = 1:1
        drawnow;
        frame = getframe(gca,FrameRegion);
        writeVideo(writerObj,frame);
    end
end

% stand
for i = 1:20
    frame = getframe(gca,FrameRegion);
    writeVideo(writerObj,frame);
end

% move northwest
viewregion_x = get(gca,'XLim');
viewregion_y = get(gca,'YLim');

for k = 1:50
    viewregion_x = viewregion_x - range(viewregion_x)*0.004;
    viewregion_y = viewregion_y + range(viewregion_y)*0.004;   
    axis([viewregion_x,viewregion_y]);
    
    for i = 1:2
        drawnow;
        frame = getframe(gca,FrameRegion);
        writeVideo(writerObj,frame);
    end
end

% stand
for i = 1:30
    frame = getframe(gca,FrameRegion);
    
    writeVideo(writerObj,frame);
end

% change to the original view (zoomed out)
zoom out;
for i = 1:20
    frame = getframe(gca,FrameRegion);
    
    writeVideo(writerObj,frame);
end

close(writerObj);