
image_directory = 'H:\Sample 1\';
files = {
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

imdir=image_directory

for s = 1:numel(files)
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    % Construct a Bio-Formats reader decorated with the Memoizer wrapper
    bfreader = loci.formats.Memoizer(bfGetReader(), 0);
    % Initialize the reader with an input file
    bfreader.setId(fullfile(indir, imfile));
    
    
end