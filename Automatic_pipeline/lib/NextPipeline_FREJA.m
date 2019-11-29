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
 


% Choose functions to run
run_MippingAlign_YN = 0;
run_Anchor_012 = 0; % 1 means generate pseudoanchor, 2 means use real anchor % 0 means skip
        anchorchannel = 2
run_BlobIdentification_YN = 1;
run_DecodingPlotting_YN = 1;

%================================================
% set parameters
%----------------------------

out_fold='F:\HB_SECTION1';
status=mkdir(out_fold);
% Stitching
me.default_folder= 'G:\HB_TEST1\Section1';
me.rawfiles={
    'Base_1.czi',
    'Base_2.czi',
    'Base_3.czi',
    'Base_4.czi'
    };

me.outputdirectory= [out_fold,'\Preprocess']
me.nChannels = 4; %Pseudoanchor is not undestood as a channel. I think is number of bases
me.reference_round = 1;
me.tile_align_channel = 1; % for registration between rounds
me.tile_stitch_channel = 1; %Cy3   % for registration between tiles
nbases = size(me.rawfiles,1) %counts the number of files introduced in rawfiles
 % this will be the number of channels
%------------------------------

% Create csv file
cs.input_directory=[out_fold,'\Preprocess\2DTiles'];
file_name='Base_1.czi';
cs.file_name = regexprep(file_name,'.czi','');
cs.output_directory=[out_fold,'\Preprocess']

%Create anchor
PATH= [me.outputdirectory,'\Aligned2DTiles\'] %stablish input directory for creating anchor as input file for anchor
OUTdir=[out_fold,'\withanchor']
OUTPATH=[OUTdir,'\']
FILE_PREFIX='190917_msSBH_02_rd6cy1(1)_t1.tif'
FILE_SUFIX=''
%----------------------------
% Spot identification

si.folder=OUTPATH;
si.PREFIX='Tile';
si.INTERFIX='_cyc_';
si.SUFIX='_CH_';
si.TYPE='.tif';
si.outputfolder=OUTdir;
si.nameoutput='blobs.csv';
si.outputdirect=[si.outputfolder,'\',si.nameoutput];

%Decoding_and_plotting
d.input_file=[out_fold,'\withanchor\blobs.csv']
d.General_Alignment_ACGT_column_number = [0,0,5,4,3,2];
d.XYPosition_ParentCell_column_number = [7,8,0];
d.num_hybs = 4;
d.taglist = 'F:\HB_SECTION2\taglist_HBTEST1.csv';   % old .m taglist or .csv file with columns: code, name, symbol(optional), no header
d.csv_file_contain_tile_position = [out_fold,'\Preprocess\tiled.csv'];
d.output_directory = [out_fold,'\Decoding'];   
q.quality_threshold = 0.7;
p.background_image = 'F:\HB_SECTION2\Preprocess\2DTiles\Base_1_t1.tif'; 
 p.I_want_to_plot_on_white_backgound = 1;
 
 

  
%----------------FROM THIS POINT NOTHING NEED S TO BE MODIFIED

  
%===================================================================
%===================================================================
%===================================================================


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
   Makeanchor_automatic_function(PATH,OUTPATH,FILE_PREFIX,FILE_SUFIX,ntiles,nbases);
   disp("Automatic anchor made")
end
if run_Anchor_012 == 2
  Noanchor_automatic_function(PATH,OUTPATH,FILE_PREFIX,FILE_SUFIX,ntiles,nbases,anchorchannel)
end

if run_BlobIdentification_YN
    CellProfilerInMatlab_function(si.folder,si.PREFIX,si.INTERFIX,si.SUFIX,si.TYPE,ntiles,si.outputdirect,nbases);
end

if run_DecodingPlotting_YN
   decoding_plotting_function(d.input_file,d.General_Alignment_ACGT_column_number,...
    d.XYPosition_ParentCell_column_number,d.num_hybs,d.taglist,d.csv_file_contain_tile_position,...
    d.output_directory,q.quality_threshold,p.background_image, p.I_want_to_plot_on_white_backgound);
end


