function decoding_plotting_function( input_file,General_Alignment_ACGT_column_number,...
    XYPosition_ParentCell_column_number,num_hybs,taglist,csv_file_contain_tile_position,...
    output_directory,quality_threshold,background_image, I_want_to_plot_on_white_backgound);

% Choose functions to run
run_Tiling_YN = 0;
run_Decode_YN = 1;
run_Threshold_YN = 1;
run_Plotting_Global_YN = 1;
%================================================
%[0,0,7,9,8,6]
    % options
    d.General_Alignment_ACGT_column_number=General_Alignment_ACGT_column_number;
    d.XYPosition_ParentCell_column_number=XYPosition_ParentCell_column_number;
    d.csv_file_contain_tile_position=csv_file_contain_tile_position;
    d.output_directory=output_directory;
    q.quality_threshold=quality_threshold;
    p.background_image=background_image;
    p.I_want_to_plot_on_white_backgound=I_want_to_plot_on_white_backgound;
    
    
    d.num_hybs=num_hybs;
    d.taglist=taglist;
    d.input_file=input_file;
    d.check_parent_cell_YN = 0;       
    d.check_alignment_YN = 0;
        d.alignment_min_threshold = 1.8;
    d.abnormal_sequencing_YN = 0;
        d.sequencing_order = '01234';  % keep the quote marks, same length as
%----------------------------
% Thresholding       
    q.general_stain_threshold = 0;
%----------------------------
        p.scale = 1;		% image scale
    % options
    p.exclude_NNNN_YN = 0;
    p.plot_reads_beforeQT_YN = 0;
    p.plot_ref_general_stain = 0; 
%================================================

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
end
