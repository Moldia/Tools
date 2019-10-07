%% Version 1 aligns all Bases with the automatically selected
% reference image for a batch of samples.
% 2018-04-20, Erik Samuelsson
% Updates: now includes a prompt that allows the user to change the
% selected refereance image if desired.
tic
    clear;
    close all;
%% Prompt user for number of channels, bases --> decides number of loops
% Prompt user for image file paths

    prompt = {'Enter the number of Bases:','Enter the number of channels:'...
        ,'Enter the file path of the folder containing all of the images for analysis:',...
        'Enter the channel number for the DAPI stain:',...
        'Enter the number of samples to be processed:','Enter the tile size',....
        'Enter Tile Channel 1','Enter Tile Channel 2','Enter Tile Channel 3',...
        'Enter Tile Channel 4','Enter Tile Channel 5','Enter Tile Channel 6'};
    dlg_title = 'Input';
    num_lines = 1;
     defaultans = {'3','6','E:\Bad Reference base selection'...
        ,'1','1','2000','Nuclei','General_stain','T','G','C','A'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    bases = str2double(cell2mat(answer(1)));
    channels = str2double(cell2mat(answer(2)));
    filepath_samples = char(answer(3)); 
    DAPI_channel = str2double(cell2mat(answer(4)));
    samples =  str2double(cell2mat(answer(5)));
    tile_size = str2double(cell2mat(answer(6)));
    tile_channel_1 = char(answer(7));
    tile_channel_2 = char(answer(8));
    tile_channel_3 = char(answer(9));
    tile_channel_4 = char(answer(10));
    tile_channel_5 = char(answer(11));
    tile_channel_6 = char(answer(12));
%%   
% E:\Testbatch3_Marcus
for s = 1:samples
    filepath_image = [filepath_samples '\' num2str(s)];
    % filepath_flo_image
%% Entropic Method of Determining Reference Image
    Ent = zeros(1,bases);
        for b = 1:bases   
            img = imread([filepath_image '\Base ' num2str(b) '_c' num2str(DAPI_channel) '_ORG.tif']);
            Ent(1,b) = entropy(img);
        end 

    [row, computedbaseRef] = find(ismember(Ent, max(Ent(:))));

% Report selected reference Base, and ask if user wants to continue
    opts.Interpreter = 'tex';
    % Include the desired Default answer
    opts.Default = 'No';
    % Use the TeX interpreter to format the question
    message = ['Base ' num2str(computedbaseRef) ' was selected to be the reference image for registration. Would you like to select a different reference Base?'];
    answer = questdlg(message,'Reference Image Selection',...
                      'Yes','No',opts);
% Handle response              
    switch answer
        case 'Yes'
            prompt = {'Which Base is your reference? (Please enter a number):'};
            dlg_title = 'Input';
            num_lines = 1;
            defaultans = {''};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            baseRef = str2double(cell2mat(answer(1)));
            file_ref_image = [filepath_image '\Base ' num2str(baseRef) '_c' num2str(DAPI_channel) '_ORG.tif'];
        case 'No'
            baseRef = computedbaseRef;
            file_ref_image = [filepath_image '\Base ' num2str(baseRef) '_c' num2str(DAPI_channel) '_ORG.tif'];
    end

% Save automatically selected Reference Image
    save([filepath_image '\Selected Reference Image.mat'],'file_ref_image');
%% Image Registration of all Bases to the selected Reference Image

R_before = zeros(bases-1,2);
R_after = zeros(bases-1,2);
for b = 1:bases
    if b == baseRef
        % make sure that the reference image does not register to itself
        % and
        % Save all channels to the Aligned Images Folder
        for c = 1:channels
            if not(exist([filepath_image,'\Aligned_Images_Rigid'],'dir'))
                mkdir([filepath_image,'\Aligned_Images_Rigid'])
            end
        flo = imread([filepath_image '\Base ' num2str(baseRef) '_c' num2str(c) '_ORG.tif']);
        imwrite(flo,[filepath_image '\Aligned_Images_Rigid\Base '...
            num2str(baseRef) '_c' num2str(c) '_ORG.tif']);
        end 
       continue
    end
%     Run script for image registration
    file_flo_image = [filepath_image '\Base ' num2str(b) '_c' num2str(DAPI_channel) '_ORG.tif'];
%% Input: Read in the image pair

    ref_image = imread(file_ref_image); % give sample snapshot image (blue DAPI)
    flo_image = imread(file_flo_image); % give sample snapshot image (blue DAPI)
    
% Extract base number from reference image
    % Assumes following format for name: Base #_c#_ORG
    
    [filepath_ref,name_ref,ext_ref] = fileparts(file_ref_image);
    ref_image_name = strsplit(name_ref,'_');
    ref_image_prefix = cell2mat(ref_image_name(1));
    ref_image_Base = strsplit(char(ref_image_name(1)));
    ref_image_Base_num = str2double(cell2mat(ref_image_Base(2)));
    
% Extract base number and channel number from floating image
    % Assumes following format for name: Base #_c#_ORG
    % Extract Base number
    
    [filepath_flo,name_flo,ext_flo] = fileparts(file_flo_image);
    flo_image_name = strsplit(name_flo,'_');
    flo_image_prefix = cell2mat(flo_image_name(1));
    output_image_prefix = flo_image_prefix;
    flo_image_Base = strsplit(char(flo_image_name(1)));
    flo_image_Base_num = str2double(cell2mat(flo_image_Base(2)));
    
% Extract channel number

    flo_image_channel = cell2mat(flo_image_name(2));
    flo_image_channel_num = str2double(flo_image_channel(2));

    %% Correlation Coefficient of the two images before registration
    % Saves Float Base # and Correlation Coefficient into a mat file R_before

    [rows_ref, cols_ref] = size(ref_image);
    float_fit = imresize(flo_image, [rows_ref cols_ref]); 
    R_before(b,1) = b;
    R_before(b,2) = corr2(float_fit,ref_image);
    
      % If the correlation coefficient of Base b to the reference is greater
    % than or equal to 0.90, there is no need to run the alignment software
    if R_before(b,2) >= 0.90 == 1
        for c = 1:channels
            if not(exist([filepath_image,'\Aligned_Images_Rigid'],'dir'))
                mkdir([filepath_image,'\Aligned_Images_Rigid'])
            end
        flo = imread([filepath_image '\Base ' num2str(b) '_c' num2str(c) '_ORG.tif']);
        imwrite(flo,[filepath_image '\Aligned_Images_Rigid\Base '...
            num2str(b) '_c' num2str(c) '_ORG.tif']);
        % Saves the R_after values as the R_before values if the the 
        % correlation coefficient of Base b to the reference is greater 
        % than or equal to 0.90
        R_after(b,1) = b;
        R_after(b,2) = R_before(b,2);
        end 
       continue
    end
    
%% Configure parameters in imregconfig
    % Since the images analyzed have different intensity distributions and are 
    % collected from multiple channels, a multimodal configuration is selected 
    %for the image registration.
    [optimizer, metric] = imregconfig('Multimodal');

    %optimizer.InitialRadius = 6.25e-3;
    optimizer.InitialRadius = 0.5e-3;
    optimizer.MaximumIterations = 200;
    metric.UseAllPixels = 0;

    R_mat = zeros(20,1);
    tformRigid_master = [];
    Rfixed = imref2d(size(ref_image));

for i = 1:20
% Registration based on Rigid(consisting of translation and rotation)

    tformRigid = imregtform(flo_image,ref_image,'similarity',optimizer,metric);

% Apply the Rigid geometric transformation output to the moving image
    % to align it with the fixed image.

    movingRegisteredRigid = imwarp(flo_image,tformRigid,'OutputView',Rfixed);
    R_after_sub = corr2(movingRegisteredRigid,ref_image);
    R_mat(i) = R_after_sub;

    if R_after_sub >= max(R_mat) == 1

        tformRigid_master = tformRigid;

    else
    end
end
    
    movingRegisteredRigid_master = imwarp(flo_image,tformRigid_master,'OutputView',Rfixed);
    % Save transformation matrix for floating base to reference base for
    % reference
    save([filepath_image '\tform_' flo_image_prefix ' registered to ' ref_image_prefix '.mat'],'tformRigid_master');
    %% Correlation Coefficient of the two images after registration
    R_after(b,1) = b;
    R_after(b,2) = corr2(movingRegisteredRigid_master,ref_image);

%% Transform images for all channels of given Base by applying the geometric
% transformation matrix to all other channels
    size_ref = size(ref_image);
    Rfixed = imref2d(size_ref);
        for c = 1:channels
            if not(exist([filepath_flo,'\Aligned_Images_Rigid'],'dir'))
            mkdir([filepath_flo,'\Aligned_Images_Rigid'])
            end
        flo = imread([filepath_ref '\' output_image_prefix '_c' num2str(c) '_ORG.tif']);
        movingRegisteredSimilarity2 = imwarp(flo,tformRigid_master,'OutputView',Rfixed);
        imwrite(movingRegisteredSimilarity2,[filepath_ref '\Aligned_Images_Rigid\'...
            output_image_prefix '_c' num2str(c) '_ORG.tif']);
        end  
  
        %% Save R_BEFORE, R_AFTER, file_ref_image,
end

    % Save R_before and R_after for reference
    save([filepath_image '\Goodness of fit before and after registration.mat'],'R_before','R_after');
 
    %% Tiling 
    % Adapted from Sequencing_v3, written by Xiaoyan, 2017

    t.folder_image = [filepath_image '\Aligned_Images_Rigid'];
    t.filename_base_prefix = 'Base ';  % keep single quote marks
    t.filename_channel_prefix = '_c';
        t.in_subfolder_YN = 0;
    t.filename_suffix = '_ORG.tif';
    t.base_start = 1;     t.base_end = bases;       
    t.channel_start = 1;  t.channel_end = channels;
    t.tile_size = tile_size;
    t.channel_order = {tile_channel_1 tile_channel_2 tile_channel_3...
        tile_channel_4 tile_channel_5 tile_channel_6};
    t.CSV_filename_prefix = 'Tiled';
    seqtiling(t);  
end