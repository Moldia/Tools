base = 5;
channel_shift = [0, 0; 
    0, 2; 
    0, 0; 
    1, -1];
channel = [6,5,4,3];
% mkdir('channel_aligned');

tic
for i = 1:4
    if channel_shift(i,1) || channel_shift(i,2)
        xleft = channel_shift(i,1);
        yup = channel_shift(i,2);
        
        for j = 1:base
            I = imread(['aligned_images\base' num2str(j) '_c' num2str(channel(i)) '_ORG.tif']);
            if yup>=0
                I = [I(yup+1:end,:);zeros(yup,size(I,2),'uint16')];
            else
                I = [zeros(-yup,size(I,2),'uint16');I(1:end+yup,:)];
            end
            if xleft>=0
                I = [I(:,xleft+1:end),zeros(size(I,1),xleft,'uint16')];
            else
                I = [zeros(size(I,1),-xleft,'uint16'),I(:,1:end+xleft)];
            end
            
            imwrite(I,['aligned_images\base' num2str(j) '_c' num2str(channel(i)) '_ORG.tif'],...
                'tiff','compression','none');
        end
    end
end
toc