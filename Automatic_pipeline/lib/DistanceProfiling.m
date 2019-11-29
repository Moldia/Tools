
figure;
for i = 1:4
    I = imread(['image' num2str(i) '.jpg']);
    I = I(:,:,3);
    I = max(I(:))-I;

    Iconv = I;
    subplot(2,2,i);
    
    imshow(Iconv);
    hold on;
    
    for z = 1:10
        Iconv = imdilate(Iconv,strel('disk',20));
        bound1 = bwboundaries(Iconv,8,'noholes');
        for i = 1:length(bound1)
            plot(bound1{i}(:,2),bound1{i}(:,1));
        end
    end
    set(gca,'YDir','reverse');
    axis image;
    title(['size: ' num2str(size(I,2)) 'x' num2str(size(I,1))])
end

for i = 1:4
    I = imread(['image' num2str(i) '.jpg']);
    I = I(:,:,3);
    I = max(I(:))-I;

    Iconv = I;
    subplot(2,2,i);
    
%     imshow(Iconv);
%     hold on;
    
    for z = 1:10
        Iconv = imerode(Iconv,strel('disk',20));
        bound1 = bwboundaries(Iconv,8,'noholes');
        for i = 1:length(bound1)
            plot(bound1{i}(:,2),bound1{i}(:,1));
        end
    end
    set(gca,'YDir','reverse');
    axis image;
    title(['size: ' num2str(size(I,2)) 'x' num2str(size(I,1))])
end