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
run_MippingAlign_YN = 1; % 1 means run mipping & align, 0 means not running it
run_Anchor_012 = 0; % 1 means generate pseudoanchor, 2 means use real anchor % 0 means skip
        anchorchannel = 5 % In case you use real anchor, specify in which channel is anchor. 
        channeldistribution=[1,2,3,4,6];%Specify,DAPI, A,C,T,G %only for Pseudoanchor
run_BlobIdentification_YN = 0; % 1 means doing blob Identification ,0 means not to run
run_DecodingPlotting_YN = 0;

%================================================
% set parameters
%----------------------------

out_fold='L:\whole_organoid_PSEUDOanchor';   %Select an output folder, that will be generated
status=mkdir(out_fold);
% Stitching
me.default_folder= 'H:\Sample 1'; % Select folder where the different Bases are 
me.rawfiles={
    'Base_1.czi',
    'Base_2.czi',
    'Base_3.czi',
    'Base_4.czi',
    'Base_5.czi'
    };

me.outputdirectory= [out_fold,'\Preprocess']
me.nChannels = 5; %Number of channels NOTE THAT:Pseudoanchor is NOT undestood as a channel. Do not count on it.
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

%Create anchor
PATH= [me.outputdirectory,'\Aligned2DTiles\'] 
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
d.input_file=[out_fold,'\withanchor\blobsNO.csv'] 
d.General_Alignment_ACGT_column_number = [0,0,3,4,2,5]; %Add columns of A ,C ,G and T in this order
d.XYPosition_ParentCell_column_number = [7,8,0]; % Add colums for X and Y position
d.num_hybs = me.nChannels;
d.taglist = 'L:\whole_organoid\Output_from_standard\taglist_Organoid_fixed_colors.csv';   % old .m taglist or .csv file with columns: code, name, symbol(optional), no header
d.csv_file_contain_tile_position = [out_fold,'\Preprocess\tiled.csv']; 
d.output_directory = [out_fold,'\Decoding'];   
q.quality_threshold = 0.7; % Modify the minimum quality
p.background_image = ['H:\iss-analysis-master\lib\omero\bfmatlab\PreprocessingERIK2\OUTPUT ALIGNMENT\Base 1_c1_ALI.tif']; %Select background image if desired
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


if run_MippingAlign_YN || run_Anchor_012 || run_BlobIdentification_YN || run_DecodingPlotting_YN
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
a=dir([out_fold,'/Preprocess/2DTiles', '/*.tif'])
elements_in_Aligned2DTiles=size(a,1)
ntiles=elements_in_Aligned2DTiles/me.nChannels


if run_Anchor_012 == 1
  Makeanchor_automatic_function(PATH,OUTPATH,FILE_PREFIX,FILE_SUFIX,ntiles,nbases,channeldistribution);
end
if run_Anchor_012 == 2
  Noanchor_automatic_function(PATH,OUTPATH,FILE_PREFIX,FILE_SUFIX,ntiles,nbases,anchorchannel)
end

if run_BlobIdentification_YN
    CellProfilerInMatlab_function(si.folder,si.PREFIX,si.INTERFIX,si.SUFIX,si.TYPE,ntiles,si.outputdirect);
end

if run_DecodingPlotting_YN
    exclude_overlapping_regions (o.pathtotiled,o.tiledcsv,o.pathtoblobs,o.blobscsv);
    decoding_plotting_function(d.input_file,d.General_Alignment_ACGT_column_number,...
    d.XYPosition_ParentCell_column_number,d.num_hybs,d.taglist,d.csv_file_contain_tile_position,...
    d.output_directory,q.quality_threshold,p.background_image, p.I_want_to_plot_on_white_backgound);
end


