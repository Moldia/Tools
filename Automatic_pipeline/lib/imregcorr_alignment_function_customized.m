
%% Version 2 (beta) aligns all Bases with the automatically selected
% 2019-09-20, Erik Samuelsson & Sergio Marco
% Developed at Molecular Diagnostic group, SciLifeLab, Stockholm University
% Updates: now includes a prompt that allows the user to change the
% selected refereance image if desired.


function imregcorr_alignment_function_customized(nbases,nchannels, input_folder,dapi,tile_size,ch1,ch2,ch3,ch4,ch5,ch6,...
    subsection,ref_base,outfolder)
    number_bases = nbases;
    number_channels = nchannels;
    filepath_samples = input_folder;
    DAPI_channel = dapi;
    samples = 1;
    tile_size = tile_size;
    tile_channel_1 =ch1;
    tile_channel_2 =ch2;
    tile_channel_3 =ch3;
    tile_channel_4 =ch4;
    tile_channel_5 =ch5;
    tile_channel_6 =ch6;
    totaltiles = subsection;
    reference_base = ref_base;
    OUTPUT_FOLDER = outfolder;     
    
%Other parameters that should't be modified   
%The gap between comparisons. gap=1 is the most precise. gap>1 gives faster approaches.
INPUT_FOLDER_REF=filepath_samples;
INPUT_FOLDER_FLO=filepath_samples;
%FLOAT_FILE='Base 1_c1.tif';

%Max pixel/gap has to be an entire number
movx=0;
movy=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Parameters present before in function%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir(OUTPUT_FOLDER)
% downscaling the images due to RAM limitations
ref_image_prefix=['Base ' num2str(reference_base) '_c'];
ref_image = imread([INPUT_FOLDER_FLO ref_image_prefix num2str(DAPI_channel) '_ORG.tif']);
ref_image_re = ref_image;

% This just saves reference base in the new format

for c = 1:number_channels
    output_image_prefix_ref=['Base ' num2str(reference_base) '_c']
    flo = imread([INPUT_FOLDER_REF ref_image_prefix num2str(c) '_ORG.tif']);
    imwrite(flo,[OUTPUT_FOLDER '\' output_image_prefix_ref num2str(c) '.tif'])

end



chan=[1:number_bases];
chan=chan(chan~=(reference_base))

COR=zeros(number_bases-1,number_channels)
for si=1:size(chan,2)
w=chan(si)    
input_image_prefix=['Base ' num2str(w) '_c']
output_image_prefix=['Base ' num2str(w) '_c']
flo_image = imread([INPUT_FOLDER_FLO input_image_prefix num2str(DAPI_channel) '_ORG.tif']);
float_fit= flo_image;
% RAM management
clear ref_image;
clear flo_image;

% resizing the floating image to be the same dimensions as the reference
% image


%Image resizing
[rows_ref, cols_ref] = size(ref_image_re);
[rows_flo, cols_flo] = size(float_fit);

tilewidth=round(0.1*((rows_ref+cols_ref)/2),0)
xtilewidth=round(0.1*((cols_ref)),0)
ytilewidth=round(0.1*((rows_ref)),0)

% YOU SHOULD GET FLOAT_FIT & REF_IMAGE_RE
% RAM management 
clear flo_image_re; 


%% 


[Xlim,Ylim]=size(float_fit);

xleftotal=[]
yuptotal=[]

tt=totaltiles;

TOTAL={}

%%
%HERE RUN AGAIN

for i =1:(tt);
for s=1:(tt);
 
[Xlim,Ylim]=size(float_fit);
xMax_width=ytilewidth;
yMax_width=xtilewidth;
ref_tile=ref_image_re(round(Xlim/(tt+1),0)*i:(round(Xlim/(tt+1),0)*i)+xMax_width,round(Ylim/(tt+1),0)*s:round((Ylim/(tt+1)),0)*s+yMax_width);
float_tile=float_fit(round(Xlim/(tt+1),0)*i:(round(Xlim/(tt+1),0)*i)+xMax_width,round(Ylim/(tt+1),0)*s:round((Ylim/(tt+1)),0)*s+yMax_width);
size(ref_tile)
size(float_tile)
if (sum(sum(float_tile))/tilewidth> 200000) & (sum(sum(ref_tile))/tilewidth> 200000)
disp(i);
tform = imregcorr(ref_tile,float_tile);   
disp(s);
%[xleft,yup] = region(ref_tile,float_tile,MAXIM_PIXEL,number_channels,input_image_prefix,OUTPUT_FOLDER,gap,movx,movy);
TOTAL{end+1}=tform;  

end 
end
end

tform=TOTAL{1};
bigm=zeros(3,3,size(TOTAL,2));
for i = 1:size(TOTAL,2)
    bigm(:,:,i)=TOTAL{i}.T;
end

for a=size(tform.T,1)
    for b=size(tform.T,2)
        tform.T(a,b)=mode(bigm(a,b,:));
    end
end


for c = 1:number_channels

    flo = imread([INPUT_FOLDER_REF input_image_prefix num2str(c) '_ORG.tif']);
    Rfixed = imref2d(size(ref_image_re));
    movingReg = imwarp(flo,tform,'OutputView',Rfixed);
    corafter=corr2(ref_image_re,movingReg);
    COR(si,c)=corafter
    imwrite(movingReg,[OUTPUT_FOLDER '\' output_image_prefix num2str(c) '.tif'])

end


end
COR=array2table(COR)
writetable(COR,[OUTPUT_FOLDER '\cross_correlations_after.csv'])







%End of the part to modify in the script




end
