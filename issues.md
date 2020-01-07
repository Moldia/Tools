issues in image analysis

- [1. javabridge issue on cp](#1-javabridge-issue-on-cp)
  * [1.1 description](#11-description)
  * [1.3 solution](#13-solution)
- [2. issues with mipping and stitching](#2-issues-with-mipping-and-stitching)
  * [2.1 description](#21-description)
  * [2.3 solution](#23-solution)
- [3. issues with the matlab path](#3-issues-with-the-matlab-path)
  * [3.1 description](#31-description)
  * [3.2 solution](#32-solution)
- [4. issue with the automatic pipeline and mipping](#4-issue-with-the-automatic-pipeline-and-mipping)
  * [4.1 description](#41-description)
  * [4.2 solution](#42-solution)
- [5. issues in the decoding - index exceeds the number of array elements](#5-issues-in-the-decoding---index-exceeds-the-number-of-array-elements)
  * [5.1 description](#51-description)
  * [5.2 solution](#52-solution)
- [6. issues with a lot of nnnns](#6-issues-with-a-lot-of-nnnns)
  * [6.1 description](#61-description)
  * [6.2 solution](#62-solution)
- [7. issues in the plotting, not seeing all the genes](#7-issues-in-the-plotting--not-seeing-all-the-genes)
  * [7.1 description](#71-description)
  * [7.2 solution](#72-solution)
- [8. issues in format_for_starfish](#8-issues-in-format-for-starfish)
  * [8.1 description](#81-description)
  * [8.2 solution](#82-solution)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

# 1. javabridge issue on cp
## 1.1 description 
this issue is quite commonly occuring and it occurs in cp after a set of images have been processed. this is dependent on the computer but it can arise after 500 images or up to 700 images for instance. 

## 1.3 solution 
the solution for this issue is to simply split up the processing meaning that we will process just a range of rows. this setting is present in cp, where there is a option to chose "Rows to process". 

after this, the files have to be put together again, which requires some finesse. 

# 2. issues with mipping and stitching
## 2.1 description 
the problem stems from the fact that we are running the old version of zen on the microscope computer and that means that mipping does not always work. 

## 2.3 solution 
the solution for this issue is to simply use zen 2.0 when we run into issues with the mipping and then run 2.3 on the computer. we need to have one computer with 2 and one with 2.3. 

# 3. issues with the matlab path

## 3.1 description
the issues comes up in matlab and it says that the default path is not working.  a

## 3.2 solution 
one way to work around this issue is to first recreate the issue by typing into the command window, "which pathdef" and then "restoredefaultpath". this will restore the matlab path to the default.

# 4. issue with the automatic pipeline and mipping
## 4.1 description 
a issue that I'm running into is that when i try to run the automatic pipeline, is that one large file is created as opposed to it actually tiling the images. it also seems to be a tad bit random as to what the image is from. in addition, when the second cycle is being processed, the output tif images seems to be the exact same as the first round image. when viewing the task manager, there are obviously things being read and written from the disk. could there be an issue with the fact that we have multiple files that are named the same? one thing that i did try was to uninstall some of the toolboxes to see if that could be the source to some of the issues. another source of the issue could be from the number of z stacks that we had. i tired running a sample set with 21 stacks compared to the 23 that we had in hca11. 

the next step was for me to try to fix the issue with running the already mipped images from zen. 

## 4.2 solution 
the solution seems to be to let the software run. it takes a long while to actually run it. one of the folders that serg created using the automatic pipeline was created 191128 @ 11:04:28, the first tile was added 191128 @ 12:25:00. the last file was modified 191291 @ 18:43:21. this seems to be the be the case with all of the implementations that make use of bioformats.  

this issue seems to be due to some issues with the format of the metadata from the microscopy images. 

# 5. issues in the decoding - index exceeds the number of array elements
## 5.1 description 
this issue arises in the last step of the analysis in which the error message reads:


``` matlab 
Index exceeds the number of array elements (0).
Error in decoding (line 139)
        globalpos(startblob+1:startblob+nBlobsInTile,1)=globalpos(startblob+1:startblob+nBlobsInTile,1)
        + tileadd(1);

Error in Tiling (line 72)
    decoding(d);
```
where the issue stems from is still unknown. there seems to be an issue with the addition of the x and y global coordinates. 

``` matlab 
%Adding x and y to global coordinates              
        globalpos(startblob+1:startblob+nBlobsInTile,1)=globalpos(startblob+1:startblob+nBlobsInTile,1) + tileadd(1);
        globalpos(startblob+1:startblob+nBlobsInTile,2)=globalpos(startblob+1:startblob+nBlobsInTile,2) + tileadd(2
```

this is the step at which it is occurs. this step is under the load and processing of the intensity data. 

## 5.2 solution 

the solution for this issue was that we had some deprecated code that had some lines that the other did not have. 

the for loop that did not work looked like follows: 

```matlab
for t = 1:max(seqdata(:,10))
    tiledata = seqdata(seqdata(:,10) == t,:);
    nBlobsInTile = size(tiledata,1)/nHybs;
    nCellsInTile = max(tiledata(:,9));
    
    fprintf('%6d\t\t%6d\n', t, nBlobsInTile);
    
    if nBlobsInTile
        [~, idxHyb1, idxBlob] = unique(tiledata(:,14));
        [~, order] = sort(idxBlob);
        
        % blob index 
        iBlob(startblob+1:startblob+nBlobsInTile, 1) =...
            (1:nBlobsInTile)' + startblob;
        
        % tile ID
        iBlob(startblob+1:startblob+nBlobsInTile, 2) = t;
        
        % max intensity and base-calling
        [maxint, maxidx] = max(tiledata(order,2:5), [], 2);
        maxChannel(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(maxint, nHybs, nBlobsInTile))';
        basecall(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(maxidx, nHybs, nBlobsInTile))';
        
        % original intensities
        sumChannel(startblob+1:startblob+nBlobsInTile,:) = ...
            reshape(sum(tiledata(order,2:5), 2), nHybs, nBlobsInTile)';
        temp = reshape(tiledata(order,2:5)', 4, nHybs, []);
        temp = permute(temp, [2, 1, 3]);
        originalChannel(:,:,startblob+1:startblob+nBlobsInTile) = temp;
        
        % alignment and general stain
        alignment(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(tiledata(order,5), nHybs, nBlobsInTile))';
        anchor(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(tiledata(order,4),nHybs,nBlobsInTile))';
        
        % global coordinates
        
         tileadd = tilepos(tilepos(:,1)==t,2:3);
        
        vector=tiledata(idxHyb1,12:13);   
          sizevect=size(vector);
          vect=[1:nHybs:sizevect(1)];
          post=tiledata(vect,12:13);
            globalpos(startblob+1:startblob+nBlobsInTile,:) =post;
        
       disp(tileadd)
        %Adding x and y to global coordinates              
        globalpos(startblob+1:startblob+nBlobsInTile,1)=globalpos(startblob+1:startblob+nBlobsInTile,1) + tileadd(1);
        globalpos(startblob+1:startblob+nBlobsInTile,2)=globalpos(startblob+1:startblob+nBlobsInTile,2) + tileadd(2);
    
        
        % parent cell
        iCell(startblob+1:startblob+nBlobsInTile) = ...
            tiledata(vect,9) + tiledata(vect,9)*startcell;        
        
        startblob = startblob + nBlobsInTile;
        startcell = startcell + nCellsInTile;
    end
    
end
```



the working for loop looked like follows: 

```matlab
for t = 1:max(seqdata(:,3))
    tiledata = seqdata(seqdata(:,3) == t,:);
    nBlobsInTile = size(tiledata,1)/nHybs;
    nCellsInTile = max(tiledata(:,12));
    
    fprintf('%6d\t\t%6d\n', t, nBlobsInTile);
    
    if nBlobsInTile
        [~, idxHyb1, idxBlob] = unique(tiledata(:,2));
        [~, order] = sort(idxBlob);
        
        % blob index 
        iBlob(startblob+1:startblob+nBlobsInTile, 1) =...
            (1:nBlobsInTile)' + startblob;
        
        % tile ID
        iBlob(startblob+1:startblob+nBlobsInTile, 2) = t;
        
        % max intensity and base-calling
        [maxint, maxidx] = max(tiledata(order,6:9), [], 2);
        maxChannel(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(maxint, nHybs, nBlobsInTile))';
        basecall(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(maxidx, nHybs, nBlobsInTile))';
        
        % original intensities
        sumChannel(startblob+1:startblob+nBlobsInTile,:) = ...
            reshape(sum(tiledata(order,6:9), 2), nHybs, nBlobsInTile)';
        temp = reshape(tiledata(order,6:9)', 4, nHybs, []);
        temp = permute(temp, [2, 1, 3]);
        originalChannel(:,:,startblob+1:startblob+nBlobsInTile) = temp;
        
        % alignment and general stain
        alignment(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(tiledata(order,5), nHybs, nBlobsInTile))';
        anchor(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(tiledata(order,4),nHybs,nBlobsInTile))';
        
        % global coordinates
        try
            globalpos(startblob+1:startblob+nBlobsInTile,:) =...
                bsxfun(@plus,...
                tilepos(tilepos(:,1)==t,2:3),...
                tiledata(idxHyb1,10:11));
        catch
            globalpos(startblob+1:startblob+nBlobsInTile,:) =...
                tiledata(idxHyb1,10:11);
        end
        
        % parent cell
        iCell(startblob+1:startblob+nBlobsInTile) = ...
            tiledata(idxHyb1,12) + tiledata(idxHyb1,12)*startcell;        
        
        startblob = startblob + nBlobsInTile;
        startcell = startcell + nCellsInTile;
    end
    
end
```
# 6. issues with a lot of nnnns
## 6.1 description 
the issue that im running into here is that we get a ton of nnnns. out of the a total of well over 200 000 genes, half of them are nnnns. whatever this issue seems to stem from, it has to be fixed. could there be a similar issues as to the one that erik ran into. 

there does not seem to be an issue with the alignment, since when i compared it to the output the serg got from the hca10, his looked similar to mine in the alignment step. 

![image-20200107101202425](C:\Users\chris.langseth\AppData\Roaming\Typora\typora-user-images\image-20200107101202425.png)

and the issue is not that im doing the same thing as erik was doing. 

## 6.2 solution 


# 7. issues in the plotting, not seeing all the genes 
## 7.1 description 
the issue arises when i try to plot the genes and some of the genes are not included in the output. 

## 7.2 solution

# 8. issues in format_for_starfish
## 8.1 description 
this issue comes up when i try to run the format_for_starfish in which the function get_2d_tiles does not recognize the tiles in the czi files. 

## 8.2 solution
