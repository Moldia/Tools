

t.folder_image = [filepath_samples '\Aligned_Images_Rigid'];
    t.filename_base_prefix = 'Base_';  % keep single quote marks
    t.filename_channel_prefix = '-aligned';
        t.in_subfolder_YN = 0;
    t.filename_suffix = '_ALI.tif';
    t.base_start = 1;     t.base_end = number_bases;       
    t.channel_start = 1;  t.channel_end = number_channels;
    t.tile_size = tile_size;
    t.channel_order = {tile_channel_1 tile_channel_2 tile_channel_3...
        tile_channel_4 tile_channel_5 tile_channel_6};
    t.CSV_filename_prefix = 'Tiled';
    seqtiling(t);  

else
    disp('Tiling will NOT be done, as specified by the user')
    
end