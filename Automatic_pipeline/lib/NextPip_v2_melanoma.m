% NEXTPIPELINE
% Developed by Sergio Marco & Xiaoyan Qian
% Next generation in situ sequencing analysis
% Ideally keep one file for one experiment
% Czi input files should be named like this: "Base_1.czi", "Base_2.czi"

% Four main functions: 
       % 1)Mipping & Align & Export: mipp files, align between channels and
       %            create .tif
       % 2) Anchor: Creates pseudoanchor in case it's needed
       % 3) BlobIdentification: identifies blobs and anotate maximum intensities in all the channels
       %        It is equivalent to CellProfiler in standard pipeline
       % 4) Decoding & Plotting: decodes the blobs.csv from previous step
       %        into gene expression
% All variables ending with _YN are yes or no (1 or 0) questions.
%
clear, clc; close all; drawnow;
 
%%%%%All parameters that need to be modified present comments

% Choose functions to run
run_MippingAlign_YN = 0; % 1 means run mipping & align, 0 means not running it
run_alignment_YN = 0;
run_tiling_YN=0;
run_Anchor_012_v = 0; % 1 means generate pseudoanchor, 2 means use real anchor % 0 means skip
        anchorchannel = 0 % In case you use real anchor, specify in which channel is anchor. 
        channeldistribution=[1,3,4,5,6];%Specify,DAPI, A,C,T,G %only for Pseudoanchor
run_BlobIdentification_YN = 0; % 1 means doing blob Identification ,0 means not to run
run_DecodingPlotting_YN = 1;

%================================================
% set parameters
%----------------------------

out_fold='E:\Melanoma';   %Select an output folder, that will be generated
status=mkdir(out_fold);
% Stitching
me.default_folder= 'J:\Melanoma #M394_1'; % Select folder where the different Bases are 
me.rawfiles={
    'Base_1.czi',
    'Base_2.czi',
    'Base_3.czi',
    'Base_4.czi',
    };

me.outputdirectory= [out_fold,'\Preprocess']
me.nChannels = 6; %Number of channels NOTE THAT:Pseudoanchor is NOT undestood as a channel. Do not count on it.
me.reference_round = 1; % Round to use as a reference
me.tile_align_channel = 1; % Align channel used 
me.tile_stitch_channel = 1; %Cy3   % for registration between tiles
nbases = size(me.rawfiles,1) 

%------------------------------

% Create csv file
cs.input_directory=[out_fold,'\Preprocess\2DTiles'];
file_name=me.rawfiles{1};
cs.file_name = regexprep(file_name,'.czi','');
cs.output_directory=[out_fold,'\Preprocess']

%Alignment parameters
nbases_ali=nbases
nchannels_ali= me.nChannels
input_folder_ali =[out_fold,'\Preprocess\Stitched2DTiles_MIST_Ref1\']
dapi_ali=1;
tile_size_ali=2000;
ch1='Nuclei'
ch2='General_stain'
ch3='T'
ch4='G'
ch5='C'
ch6='A'
subsection_ali=5
ref_base_ali=3
outfolder_ali=[out_fold,'\Align']


%Tiling parameters
t.folder_image = [out_fold,'\Align'];
    t.filename_base_prefix = 'Base_';  % keep single quote marks
    t.filename_channel_prefix = '_aligned-';
        t.in_subfolder_YN = 0;
    t.filename_suffix = '.tif';
    t.base_start = 1;     t.base_end = nbases;       
    t.channel_start = 1;  t.channel_end = me.nChannels;
    t.tile_size = 100; %Modify tile size
    t.channel_order = {'Nuclei', 'G', 'A', 'C', 'T','General_stain'};
    t.CSV_filename_prefix = [out_fold,'\Align\tiledpos'];
    t.in_subfolder_YN=0





%Create anchor
PATH= [out_fold,'\Align\'] 
OUTdir=[out_fold,'\withanchor']
OUTPATH=[OUTdir,'\']
FILE_PREFIX=''
FILE_SUFIX=''
%----------------------------
% Spot identification

si.folder=OUTPATH;
si.PREFIX='Tile';
si.INTERFIX='_cyc_';
si.SUFIX='_CH_';
si.TYPE='.tif';
si.outputfolder=OUTdir;
si.nameoutput='blobs.csv'; % Select the name of blobs file
si.outputdirect=[si.outputfolder,'\',si.nameoutput]; %LOCATION OF blobs (shouldn't be modified)

%--------------------------------
%Exclude overlap
o.pathtotiled = [out_fold,'\Preprocess\']
o.tiledcsv = 'tiled.csv'
o.pathtoblobs = [out_fold,'\withanchor\']
o.blobscsv = 'blobs.csv'
%--------------------------------
%Decoding_and_plotting
d.input_file=[out_fold,'\withanchor\blobs.csv'] 
d.General_Alignment_ACGT_column_number = [0,0,3,4,2,5]; %Add columns of A ,C ,G and T in this order
d.XYPosition_ParentCell_column_number = [7,8,0]; % Add colums for X and Y position
d.num_hybs = nbases;
d.taglist = [out_fold,'\taglist.csv']; % old .m taglist or .csv file with columns: code, name, symbol(optional), no header

d.csv_file_contain_tile_position=[out_fold,'\Align\tiledpos.csv']; 
d.output_directory = [out_fold,'\Decoding'];   
q.quality_threshold = 0.4; % Modify the minimum quality
p.background_image = [out_fold,'\Align\Base_1_aligned-1.tif']; %Select background image if desired
 p.I_want_to_plot_on_white_backgound = 0; % Select whether a background image is needed (1) or not (0)
 
 

  
%----------------FROM THIS POINT NOTHING NEED S TO BE MODIFIED

%=====================================================================================================
%=====================================================================================================
%=====================================================================================================


%Make all the directories related.
mkdir(me.outputdirectory)
mkdir(OUTPATH)
mkdir(si.outputfolder)
mkdir(d.output_directory)


if run_MippingAlign_YN || run_Anchor_012_v || run_BlobIdentification_YN || run_DecodingPlotting_YN ||run_tiling_YN||run_alignment_YN
else
    error('Choose at least one function.');
end
%Run czi exporting to mipping using bioformats and alignment on each tile
if run_MippingAlign_YN
  mippingexport_function(me.default_folder,me.rawfiles,me.outputdirectory,me.nChannels,me.reference_round,me.tile_align_channel,...
      me.tile_stitch_channel)
  d.csv_file_contain_tile_position = adapt_tiles(cs.input_directory,cs.file_name,cs.output_directory)
end

% Calculate 
a=dir([out_fold,'/Align/Tiled_Base_1_aligned-1', '/*.tif'])
elements_in_Aligned2DTiles=size(a,1)
ntiles=elements_in_Aligned2DTiles

if run_alignment_YN == 1;
    imregcorr_alignment_function(nbases_ali,nchannels_ali, input_folder_ali,...
        dapi_ali,tile_size_ali,ch1,ch2,ch3,ch4,ch5,ch6,...
    subsection_ali,ref_base_ali,outfolder_ali)
    
end
if run_tiling_YN == 1;
     seqtiling(t)
end


if run_Anchor_012_v == 1;
  Makeanchor_automatic_function(PATH,OUTPATH,FILE_PREFIX,FILE_SUFIX,ntiles,nbases,channeldistribution);
end
if run_Anchor_012_v == 2;
  Noanchor_automatic_function(PATH,OUTPATH,FILE_PREFIX,FILE_SUFIX,ntiles,nbases,anchorchannel)
end


if run_BlobIdentification_YN;
    CellProfilerInMatlab_function_2(si.folder,si.PREFIX,si.INTERFIX,si.SUFIX,si.TYPE,ntiles,si.outputdirect,nbases);
end

if run_DecodingPlotting_YN
    %exclude_overlapping_regions (o.pathtotiled,o.tiledcsv,o.pathtoblobs,o.blobscsv);
    decoding_plotting_function(d.input_file,d.General_Alignment_ACGT_column_number,...
    d.XYPosition_ParentCell_column_number,d.num_hybs,d.taglist,d.csv_file_contain_tile_position,...
    d.output_directory,q.quality_threshold,p.background_image, p.I_want_to_plot_on_white_backgound);
end


