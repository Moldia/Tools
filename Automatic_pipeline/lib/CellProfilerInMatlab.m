% ISS image analysis workshop, 170614
% Xiaoyan
% tested on MATLAB R2016b


folder='I:\3x3tile\Preprocess\withANCHOR\'
PREFIX='Tile'
INTERFIX='_cyc_'
SUFIX='_CH_'
TYPE='.tif'


TOTAL=[folder,PREFIX,num2str(tile),INTERFIX,num2str(cyc),SUFIX,num2str(channel),TYPE]




% % write image and prepare CP input file
% for b = 1:4
%     for c = 1:6
%         I = Imip{b,c};
%         imwrite(I, ['b' num2str(b) '_c' num2str(c) '.tif']);
%     end
% end
% writeposfile('b1_c1.tif', 600, 1:4, 1:6, [1,2], 'b', '_c', '', '');

%% sequencing images
% figure
figure(1); clf;
set(gcf, 'name', 'Sequencing images', 'units', 'normalized', 'position', [0 0 1 1]); 
channels = {'DAPI', 'anchor', 'T', 'G', 'C' 'A'};
Ax = [];
tile=1
for b = 1:4
    for c = 1:6
        TOTAL=[folder,PREFIX,num2str(tile),INTERFIX,num2str(b),SUFIX,num2str(c),TYPE]
        im(TOTAL)
        ax = subplot(4,7,(b-1)*7+c); Ax = [Ax, ax];
        imshow(, []);
        title(['base' num2str(b) ' ' channels{c}]); drawnow
    end
end
linkaxes(Ax, 'xy');
pause()

% add reference blob image
for b = 1:4
    ax = subplot(4,7,7*b); Ax = [Ax, ax];
    imshow(Imip{1,2}, []); title('reference anchor');
end
linkaxes(Ax, 'xy');



%% FFT phase correlation alignment
tform = imregcorr(Imip{1,2}, Imip{3,2});
subplot(4,7,16);
I = Imip{3,2};
I = padimg(I, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
imshow(I, []); title('aligned base3 anchor'); drawnow
pause()

Ialigned = cell(4,6);
Ialigned(1,:) = Imip(1,:);
% apply to all from base3
for c = 1:6
    I = Imip{3,c};
    I = padimg(I, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
    Ialigned{3,c} = I;
    subplot(4,7,14+c);
    imshow(I, []);
    title(['aligned base3 ' channels{c}]); drawnow
end
pause()

% base2 and base4
for b = 2:4
    if b ~= 3
        tform = imregcorr(Imip{1,2}, Imip{b,2});
        
        for c = 1:6
            I = Imip{b,c};
            I = padimg(I, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
            Ialigned{b,c} = I;

            subplot(4,7,(b-1)*7+c);
            imshow(I, []);
            title(['aligned base' num2str(b) ' ' channels{c}]);
        end
    end
end
drawnow
save WS170614 Ialigned -append

%% tophat of aligned
Itop = cell(4,6);
for b = 1:4 
    for c = 1:6
        I = Ialigned{b,c};
        if c == 1
            I = imtophat(I, strel('disk', 15));
        else
            I = imtophat(I, strel('disk', 2));
        end
        Itop{b,c} = I;
        subplot(4,7,(b-1)*7+c);
        imshow(I, []);
        title(['base' num2str(b) ' ' channels{c}]);
    end
end
for b = 1:4
    subplot(4,7,7*b);
    imshow(Itop{1,2}, []); title('reference anchor');
end
drawnow
save WS170614 Itop -append

%% reference segmentation
I = double(Itop{1,2})/65535;
I = imtophat(I, strel('disk', 2));
Ibw = im2bw(I, .0025);
Iws = watershed(-I);
Iws = double(Iws) .* double(Ibw);
Ibw = logical(Iws);
Ibw = bwareaopen(Ibw, 3);
Ibound = imdilate(Ibw, strel('disk', 1)) ~= Ibw;
Ibound = cat(3, uint8(Ibound*255), zeros(size(Ibound), 'uint8'), zeros(size(Ibound), 'uint8'));

subplot(4,7,7);
addimlayer(Ibound, .5);
title('reference anchor segemented');


