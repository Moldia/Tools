
% in situ sequencing pipeline
% Ideally keep one file for one experiment
% Four main functions: 1)Tiling, 2)Decoding, 3)Thresholding, 4)Plotting 
% All variables ending with _YN are yes or no (1 or 0) questions.
%
% Sequencing v3
% Xiaoyan, 2018
clear, clc; close all; drawnow;
 


% Choose functions to run
run_Tiling_YN = 0;
run_Decode_YN = 1;
run_Threshold_YN = 1;
run_Plotting_Global_YN = 1;
%================================================
% set parameters
%----------------------------
% Tiling

%----------------------------
% Decoding
%[0,0,7,9,8,6]
    d.input_file = 'L:\3x3tile\Preprocess\withANCHOR\mytable7s.csv';
    d.General_Alignment_ACGT_column_number = [0,0,4,5,2,3];    % use 0 if any of them is MISSING in the file
    d.XYPosition_ParentCell_column_number = [7,8,0];
    
    d.num_hybs = 5;
    d.taglist = 'J:\MIPs_ms120\CellProfiler\taglist_mouse120.csv';   % old .m taglist or .csv file with columns: code, name, symbol(optional), no header
    d.csv_file_contain_tile_position = 'L:\3x3tile\Preprocess\2DTiles\tiled.csv';
    d.output_directory = 'L:\3x3tile\Preprocess\Decoding';   
    % options
    d.check_parent_cell_YN = 0;       
    d.check_alignment_YN = 0;
        d.alignment_min_threshold = 1.8;
    d.abnormal_sequencing_YN = 0;
        d.sequencing_order = '01234';  % keep the quote marks, same length as
%----------------------------
% Thresholding
    q.quality_threshold = 0.7;        
    q.general_stain_threshold = 0;
%----------------------------
% Plotting
    p.background_image = 'L:\3x3tile\General pipeline process\190913_msSBH_02_rd4cy4(1)-mip_c1_ORG.tif'; 
        p.scale = 1;		% image scale
    p.I_want_to_plot_on_white_backgound = 0;
    % options
    p.exclude_NNNN_YN = 0;
    p.plot_reads_beforeQT_YN = 0;
    p.plot_ref_general_stain = 0; 
%================================================

%dirty things


%end of dirt

if run_Tiling_YN || run_Decode_YN || run_Threshold_YN || run_Plotting_Global_YN
else
    error('Choose at least one function.');
end

if run_Decode_YN
    decoding(d);
end
if run_Threshold_YN
    qthreshold(d.output_directory, q);
end 
if run_Plotting_Global_YN
    seqplotting(d.output_directory, d.taglist, q, p);
end

