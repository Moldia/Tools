
function CellProfilerInMatlab_function(folder,PREFIX,INTERFIX,SUFIX,TYPE,tiles,outputdirect)

% ISS image analysis workshop, 170614
% Xiaoyan
% tested on MATLAB R2016b




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
for tile = 1:tiles
tile
A={};
ATOT={};
%AUTOMATE
for b = 1:4
    for c = 1:6
        TOTAL=[folder,PREFIX,num2str(tile),INTERFIX,num2str(b),SUFIX,num2str(c),'_t',num2str(tile),TYPE];
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

%% FFT phase correlation alignment
%tform = imregcorr(ATOT(:,:,1,1), ATOT(:,:,1,3));
%subplot(4,7,16);
I = ATOT(:,:,6,3);
%%imshow(I, []); title('aligned base3 anchor'); drawnow

Ialigned = cell(5,6);

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
 %   I = padimg(I, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
    Ialigned{3,c} = I;
 %   subplot(4,7,14+c);
   % imshow(I, []);
    %title(['aligned base3 ' channels{c}]); drawnow
end



% base2 to base 5
for b = 2:4
    if b ~= 3
  %      tform = imregcorr(ATOT(:,:,1,1), ATOT(:,:,1,b));
        
        for c = 1:6
            I = ATOT(:,:,c,b);
  %          I = padimg(I, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
            Ialigned{b,c} = I;

            %subplot(4,7,(b-1)*7+c);
            %imshow(I, []);
            %title(['aligned base' num2str(b) ' ' channels{c}]);
        end
    end
end




%% tophat of aligned
Itop = cell(5,6);
for b = 1:4 
    for c = 1:6
        I = Ialigned{b,c};
        if c == 1
            I = imtophat(I, strel('disk', 15));
        else
            I = imtophat(I, strel('disk', 2));
        end
        Itop{b,c} = I;
        %subplot(4,7,(b-1)*7+c);
        %imshow(I, []);
        %title(['base' num2str(b) ' ' channels{c}]);
    end
end
for b = 1:4
    %subplot(4,7,7*b);
    %imshow(Itop{1,2}, []); title('reference anchor');
end


% %% reference segmentation
I = double(Itop{1,6})/65535;
I = imtophat(I, strel('disk', 2));
Ibw = im2bw(I, .0025);
Iws = watershed(-I);
Iws = double(Iws) .* double(Ibw);
Ibw = logical(Iws);
Ibw = bwareaopen(Ibw, 3);
image(Iws);
Ibw = logical(Iws);

UNIQUES = bwconncomp(Iws);
U2=labelmatrix(UNIQUES);
 blobs = unique(U2);
 blobs = blobs(blobs~=0);
 size(blobs);
 image(Iws);
%Image normalization
 for chan=1:6
    for cyc=1:4
    Itop{cyc,chan}=mat2graymod(Itop{cyc,chan});
    end
 end
 
 for totblob = 1:size(blobs)
 [X,Y]=find(U2==blobs(totblob));
 valuetottot=[];
    basetottot=[];
    channeltottot=[];
     Xpostot=[];
     Ypostot=[];
     itot=[];
     tilestot=[];
    for bases = 1:4
        valuetot=[];
        basetot=[];
        channelstot=[];
        for channels = 2:6
            fig=Itop{bases,channels};
            int=[];
            for pixels = 1:size(X,1);
                if X(pixels)<size(fig,1) & Y(pixels)<size(fig,2)
                int=[int,fig(X(pixels),Y(pixels))];  
                end
            end
            if max(int)>0
            value=max(int);
            else 
            value=0;
            end 
            valuetot = [valuetot,value];
            basetot=[basetot,bases];
            channelstot=[channelstot,channels];
            
            
            %end
        end
     basetottot=[basetottot,bases];
     tilestot=[tilestot,tile];
     valuetottot=[valuetottot,valuetot']; 
     channeltottot=[channeltottot,channelstot'];
     Xpostot=[Xpostot,mean(X)];
     Ypostot=[Ypostot,mean(Y)];
     itot=[itot,i];
 
    
    end
    
    combinat=[basetottot',valuetottot',Xpostot',Ypostot',itot',tilestot'];
    comb=size(combinat);
    if comb(2)==10
    combinattot=vertcat(combinattot,combinat);
    end
 end
end
 
CT=array2table(combinattot);

CT.Properties.VariableNames = {'ImageNumber' 'Intensity_MaxIntensity_aligned_A' 'Intensity_MaxIntensity_aligned_C' ...
    'Intensity_MaxIntensity_aligned_G' 'Intensity_MaxIntensity_aligned_T' 'Intensity_MaxIntensity_Aligned_SpecBlob' 'Location_Center_X' 'Location_Center_Y' 'ObjectNumber' 'Metadata_position'};


sizeCT=size(CT);
emty=zeros(sizeCT(1),1);
emty(1:sizeCT)=2;
CT.Intensity_MaxIntensity_AlignmentMath=emty;


writetable(CT,outputdirect)

 
 %  imshow(Ibo)
%  Ibound = imdilate(Ibw, strel('disk', 1)) ~= Ibw;
%  Ibo = bwconncomp(Ibound);
%  Ibound2 = cat(3, uint8(Ibound*255), zeros(size(Ibound), 'uint8'), zeros(size(Ibound), 'uint8'));
%  Ibound3=cat(3, uint8(Ibound*255), detection2,zeros(size(Ibound), 'uint8'));
%  subplot(4,7,7)
%  addimlayer(Ibound3, .5);
%  title('reference anchor segemented');
%  



end