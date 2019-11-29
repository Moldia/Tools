% stitch aligned images
for b = 1:3
    im = stitch_buffered_tiles...
        ('..\Preprocessing\Stitched\Ch030_170310_Pw_5dpa_SBL2_CX_c1_stitched.tif',...
        'Alignment\blobs_align', 2000, 100, '.png', 2, b:3:60);
    imwrite(im, ['Alignment_b' num2str(b) '.tif']);
end

% process blobs files
remove_overlap_blobs('blobs.csv', 2000, 100);


% % crop tile label images
% names = {'DAPI', 'Blob', 'Kv21'};
% for b = 1:3
%     outdir = ['..\Tiled_Ab_Aligned_c', num2str(b), '_ORG'];
%     mkdir(outdir);
%     for t = 1:204
%         im = imread(['Aligned\', names{b}, '_', paddigits(t,3), '.tif']);
%         im = im(51:2050, 51:2050);
%         imwrite(im,...
%             [outdir, '\tile', num2str(t), '.tif'])
%     end
% end

