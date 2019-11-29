%% MIST stitching
function stitch_2d_mist(indir, files, outdir, nChannels, ref_channel, ref_file)

outdir_top = outdir;
if nargin < 6 || ref_file==0
    ref_file = 0;
    outdir = fullfile(outdir, 'Stitched2DTiles');
else
    outdir = fullfile(outdir, ['Stitched2DTiles']);
end

mkdir(outdir);

s_order = 1:numel(files);
if ref_file
    s_order = [ref_file, setdiff(s_order, ref_file)];
    
    % transformation from tile alignment
    load(fullfile(outdir_top, 'Aligned2DTiles', 'transformation.mat'));
    TForm = cat(2, TForm{:});
    TForm = reshape(TForm, size(TForm,1), 2, []);
    aligntform = zeros(length(TForm), numel(files), 2);
    aligntform(:,:,1) = squeeze(TForm(:,1,:))';
    aligntform(:,:,2) = squeeze(TForm(:,2,:))';
    % reference image
    aligntform(:,:,1) = bsxfun(@minus, aligntform(:,:,1), aligntform(:,ref_file,1)); 
    aligntform(:,:,2) = bsxfun(@minus, aligntform(:,:,2), aligntform(:,ref_file,2)); 
    
end

save_stitched_image = 1;

imwrite(zeros(2048, 'uint16'), [indir '\empty.tif']);

for s = s_order
    
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    
    xypos = csvread(fullfile(outdir_top, '2DTiles', ['tile_coordinates_' output_prefix '.csv']));
    
    [gridorder, tilefiles] = tilepos2grid(xypos,...
        [output_prefix '_t'], '.tif', 0);
    
    output_prefix = [output_prefix '_'];

    if ref_file && s~=ref_file
        % use exisiting tramsformation
        imfile = files{ref_file};
        [~, output_prefix_ref] = fileparts(imfile);
        output_prefix_ref = [output_prefix_ref '_'];
        
        assemble_from_metadata = 1;     % apply the same to other images
        ref_transform = fullfile(outdir, [output_prefix_ref 'metadata-' num2str(ref_channel) '.mat']);
        
        for c = [ref_channel, setdiff(1:nChannels, ref_channel)]
            copyfile(ref_transform, fullfile(outdir, [output_prefix 'metadata-' num2str(ref_channel) '.mat']));
            
            % modify transformation based on tile alignment
            load(ref_transform, 'global_x_img_pos', 'global_y_img_pos');
            [tilenum, idx] = sort(gridorder(:));
            
            global_x_img_pos(idx(tilenum~=0)) = global_x_img_pos(idx(tilenum~=0)) + aligntform(:,s,1);
            global_y_img_pos(idx(tilenum~=0)) = global_y_img_pos(idx(tilenum~=0)) + aligntform(:,s,2);
            save(fullfile(outdir, [output_prefix 'metadata-' num2str(ref_channel) '.mat']),...
                'global_x_img_pos', 'global_y_img_pos', '-append');
            
            stitch_time_slice(indir, tilefiles,...
                outdir, output_prefix,...
                ref_channel, NaN, NaN, 'Overlay', 0,... 
                save_stitched_image, assemble_from_metadata,...
                fullfile(outdir, ['log_' output_prefix num2str(c) '.txt']),...
                10, 10, c);
            
            % crop/pad upper left part because starting position is always
            % set to (1,1) when assembling images for convenience
            try
                I = imread(fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif']));
                I = padimg(I, round(min(global_x_img_pos(:)))-1, round(min(global_y_img_pos(:)))-1, 'NW');
                imwrite(I, fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif']));
            catch ME
                if contains(ME.identifier, 'fileDoesNotExist')
                    load(fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif.mat']));
                    I = padimg(I, round(min(global_x_img_pos(:)))-1, round(min(global_y_img_pos(:)))-1, 'NW');
                    save(fullfile(outdir, [output_prefix 'stitched-' num2str(c) '.tif.mat']), 'I', '-append');
                end
            end
            clear I

        end
        
    else
        for c = [ref_channel, setdiff(1:nChannels, ref_channel)]
            if c == ref_channel
                % always restitch reference round reference channel
                assemble_from_metadata = 0;
            else
                assemble_from_metadata = 1;
            end
            
            stitch_time_slice(indir, tilefiles,...
                outdir, output_prefix,...
                ref_channel, 7, 2, 'Overlay', 0,...
                save_stitched_image, assemble_from_metadata,...
                fullfile(outdir, ['log_' output_prefix num2str(c) '.txt']),...
                10, 10, c);
            
        end
    end
end
end
