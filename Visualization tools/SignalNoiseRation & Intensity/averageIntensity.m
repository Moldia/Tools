% AVERAGE INTENSITY DATA 

% intensity data generated using Xiaoyans scipt called 
% PrototypeSNRProfiler.m

% organize the files as follows:
% one folder for each roundcycle. in the folder, for our purposes 
% we will have 6 subfolders, one for each ROI (X3) and sequcing chemsitry (2X)
% the naming convention for the subfodlers: 
% SBH_Round1Cycle1_SBH_ROI1
% the naming convention of the files: 
% SBH_Round1Cycle1_SBH_ROI1_c2_ORG\

clear
close
roundNumber = 6; % enter the round number 
cycleNumber = 3; % enter the cycle number

for  s = 2:5 % four channels 
    for k = ["SBH","SBL"] % two sequencing chemistries 
       Mean = [];
       Std = [];
       Tabledata = [];
       Mean_SBH = [];
       Mean_SBL = [];
       Std_SBH = []; 
       Std_SBL = []; 
        for r = ["ROI1", "ROI2", "ROI3"] % three regions of interest
                tabledata = readtable(strcat('F:\imagesAndCSVFilesFor3ROI\', 'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber), '\', (k), '_', 'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber), '_', (k), '_', (r), '\', (k), '_', 'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber), '_', (k), '_', (r), '_c', num2str(s), '_ORG.csv'));
                
%                 REMOVEDOUTLIERS = rmoutliers(tabledata, 'grubbs');
%                 REMOVEDOUTLIERS = table2array(REMOVEDOUTLIERS);
%                 IntensityDataMean = mean(REMOVEDOUTLIERS, 1); 
%                 IntensityDataStd = std(REMOVEDOUTLIERS, 1);
                
                %IntensityData = cat(1, IntensityData, middlePixelIntensity);
%                 IntensityDataMean = cat(1, Mean, IntensityDataMean);
% 
%                 
%                 IntensityDataStd = cat(1, Std, IntensityDataStd);
                
                if r == "ROI1"
                        tabledataROI1 = tabledata;
                elseif r == "ROI2" 
                        tabledataROI2 = tabledata;
                else 
                        tabledataROI3 = tabledata;
                end
 
        end
        
        if k == "SBH" 
                allDataSBH = [tabledataROI1; tabledataROI2; tabledataROI3];
                allDataSBH_rmoutliers = rmoutliers(allDataSBH); 
                meanSBH = mean(table2array(allDataSBH_rmoutliers), 1);
                stdSBH = std(table2array(allDataSBH_rmoutliers), 1);
        else 
                allDataSBL = [tabledataROI1; tabledataROI2; tabledataROI3];
                allDataSBL_rmoutliers = rmoutliers(allDataSBL);
                meanSBL = mean(table2array(allDataSBL_rmoutliers), 1);
                stdSBL = std(table2array(allDataSBL_rmoutliers), 1);
 
        end
    
    end 
    
    if s == 2
        allSBH_c2 = allDataSBH_rmoutliers;
        allSBL_c2 = allDataSBL_rmoutliers;
    elseif s == 3
        allSBH_c3 = allDataSBH_rmoutliers;
        allSBL_c3 = allDataSBL_rmoutliers;
    elseif s == 4
        allSBH_c4 = allDataSBH_rmoutliers;
        allSBL_c4 = allDataSBL_rmoutliers;
    else
        allSBH_c5 = allDataSBH_rmoutliers;
        allSBL_c5 = allDataSBL_rmoutliers;
    end    
    co = [0.0 0.0 1.0;...
          1.0 0.0 0.0];
    set(groot,'defaultAxesColorOrder',co)
    
    [h,p,ci,stats] = ttest2(meanSBH,meanSBL);
 
    sgt = sgtitle({['Average intensity over RCPs: ',  'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber)]; ''});
    sgt.FontSize = 50;
    subplot(2,2, s-1);
    errorbar(meanSBH, stdSBH, 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold on 
    errorbar(meanSBL, stdSBL, 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')
    
    tit = title(['Channel  ', num2str(s)]);
    tit.FontSize = 20 ;
    ylab = ylabel('Intensity');
    ylab.FontSize = 20; 
    
    xlab = xlabel({'Pixel number'; strcat('p value:', '    ', num2str(round(p, 2)))});
    xlab.FontSize = 20;
    
    ylim([0,    6000]);
    xlim([1,21]);
    
    lgd = legend ('SBH', 'SBL');
    lgd.FontSize = 15;         


    
end
