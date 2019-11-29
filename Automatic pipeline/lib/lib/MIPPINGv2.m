
image_directory = 'H:\Sample 1\';
image_files = {
    '190627_org1_cyc1(1).czi',
    '190629_org1_cyc2(1).czi',
    '190701_org1_cyc3(1).czi',
    '190702_org1_cyc4(1).czi',
    '190715_org1_cyc5(1).czi',
    };
output_directory = 'PreprocessingERIKOrg_v2';
nChannels = 6;
reference_round = 1;
tile_align_channel = 1; % for registration between rounds
tile_stitch_channel = 3; %Cy3   % for registration between tiles

% extract tile images from Zeiss files
%%[~, nFiles] = get_2d_tiles(image_directory, image_files, output_directory);

% align between rounds (assuming position offset is relatively small)
%%align_2d_tiles(fullfile(output_directory, '2DTiles'), image_files, output_directory, nChannels, max(nFiles), reference_round, tile_align_channel);

% stitch into big images either individually for every round or compute
% from a fixed round (need tile transformation info from align_2d_tiles)
%%stitch_2d_mist(fullfile(output_directory, '2DTiles'), image_files, output_directory, nChannels, tile_stitch_channel, reference_round);

stitch_2d_iss(fullfile(output_directory, '2DTiles'), image_files, output_directory, nChannels, tile_stitch_channel, reference_round);

% reslice into non-overlapping small tiles
nTiles = tile_stitched_images(fullfile(output_directory, 'Stitched2DTiles_MIST_Ref1'), image_files, output_directory, 2000);









