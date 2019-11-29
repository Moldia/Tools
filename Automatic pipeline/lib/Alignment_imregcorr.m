% ISS image analysis workshop, 170614
% Xiaoyan
% tested on MATLAB R2016b


folder='J:\whole_organoid_PSEUDOanchor\Preprocess\Stitched2DTiles_MIST_Ref1\'
PREFIX='Base_'
INTERFIX='_stitched-'
SUFIX=''
TYPE='.tif'

tiles = 1


%%%%%%%%%%%%%%%%%THIS VISUALIZES INFORMATION%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%CREATE A BIG OBJECT AND WORKING%%%%%%%%%%%%%%%%%%%%%%%%
combinattot=[];
A={};
ATOT={};
%AUTOMATE


for b = 1:5
    for c = 1:5
        TOTAL=[folder,PREFIX,num2str(b),INTERFIX,SUFIX,num2str(c),TYPE];
        tot=imread(TOTAL);
        if c ~= 1
            A=cat(3,A,tot);
        else 
            A=tot;
        end
       
    end
    if b ~= 1
       ATOT=cat(4,ATOT,A);
        else 
            ATOT=A;
    end
end
%last position is cycle, 3rd one is channel, previous one is y and fist is
%x

%%BLOBS FINDING

 
corbefore=corr2(tot11, tot21);
disp(corbefore);
tform = imregcorr(tot11, tot21);
I = tot2;
Rfixed = imref2d(size(tot1));
movingReg = imwarp(I,tform,'OutputView',Rfixed);
imshowpair(tot1,movingReg,'falsecolor');
corafter=corr2(tot1,movingReg);

%Base 1 saved
for c = 1:6
    I = ATOT(:,:,c,1);
    Ialigned{1,c} = I;
    %subplot(4,7,14+c);
    %imshow(I, []);
    %title(['aligned base1 ' channels{c}]); drawnow
end

% apply to all from base3
for c = 1:6
    I = ATOT(:,:,c,3);
    I = padimg(I, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
    Ialigned{3,c} = I;
    subplot(4,7,14+c);
    imshow(I, []);
    %title(['aligned base3 ' channels{c}]); drawnow
end



% base2 and base4
for b = 2:5
    if b ~= 3
        tform = imregcorr(ATOT(:,:,6,1), ATOT(:,:,6,b));
        
        for c = 1:6
            I = ATOT(:,:,c,b);
            I = padimg(I, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
            Ialigned{b,c} = I;

            %subplot(4,7,(b-1)*7+c);
            %imshow(I, []);
            %title(['aligned base' num2str(b) ' ' channels{c}]);
        end
    end
end




writetable(CT,'L:\3x3tile\Preprocess\withANCHOR\mytable7s.csv')

 
 %  imshow(Ibo)
%  Ibound = imdilate(Ibw, strel('disk', 1)) ~= Ibw;
%  Ibo = bwconncomp(Ibound);
%  Ibound2 = cat(3, uint8(Ibound*255), zeros(size(Ibound), 'uint8'), zeros(size(Ibound), 'uint8'));
%  Ibound3=cat(3, uint8(Ibound*255), detection2,zeros(size(Ibound), 'uint8'));
%  subplot(4,7,7)
%  addimlayer(Ibound3, .5);
%  title('reference anchor segemented');
%  