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
        for r = ["ROI1", "ROI2", "ROI3"] % three regions of interest
                tabledata = readtable(strcat('/Volumes/SBHvsSBL/imagesAndCSVFilesFor3ROI/', 'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber), '/', (k), '_', 'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber), '_', (k), '_', (r), '/', (k), '_', 'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber), '_', (k), '_', (r), '_c', num2str(s), '_ORG.csv'));

                REMOVEDOUTLIERS = rmoutliers(tabledata, 'grubbs');
                REMOVEDOUTLIERS = table2array(REMOVEDOUTLIERS);
                IntensityDataMean = mean(REMOVEDOUTLIERS, 1); 
                IntensityDataStd = std(REMOVEDOUTLIERS, 1);
                
                %IntensityData = cat(1, IntensityData, middlePixelIntensity);
                IntensityDataMean = cat(1, Mean, IntensityDataMean);

                
                IntensityDataStd = cat(1, Std, IntensityDataStd);

        end
        
     
        if k == "SBH" 
            % intensity mean values
            IntensityDataMeanSBH = IntensityDataMean; 
            
            % intensity standard deviation values
            IntensityDataStdSBH = IntensityDataStd;

        else 
            % intensity mean values
            IntensityDataMeanSBL = IntensityDataMean; 
            
            % intensity standard deviation values
            IntensityDataStdSBL = IntensityDataStd;

        end
    
    end 
    
    
    co = [0.0 0.0 1.0;...
          1.0 0.0 0.0];
    set(groot,'defaultAxesColorOrder',co)
   
    sgt = sgtitle(['Average intensity ',  'Round', num2str(roundNumber), 'Cycle', num2str(cycleNumber)]);
    sgt.FontSize = 30;
    subplot(2,2, s-1);
    errorbar(IntensityDataMeanSBH, IntensityDataStdSBH, 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold on 
    errorbar(IntensityDataMeanSBL, IntensityDataStdSBL, 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')
    tit = title(['Channel  ', num2str(s)]);
    tit.FontSize = 25;
    ylab = ylabel('Intensity');
    ylab.FontSize = 20; 
    ylim([0, 5000]);
    
end