% PreprocessLJStep2.m
% This file reads allData.mat file. Detrending of low frequency photometry
% signal

close all
clc  % Clear the MATLAB command window
clear all % Clear the MATLAB variables

% set parameters
params.samplerate=2052; % Hz
params.spectWindow=200; % window size for frequency calculation in spectrogram (number of samples) =97.5ms
params.spectOverlap=180; % overlap between windows in spectrogram; new calculation every 20 samples.data points thus downsampled to 102.6Hz.
params.freqRange1 = [163.5:5:178.5]; % frequencies for channel 1 spectrogram (target 171 Hz) 
params.freq1=171;
params.freqRange2 = [220.5:5:235.5]; % frequencies for channel 2 spectrogram (target 228 Hz)
params.freq2=228;
params.useFreqRange=[0:5:500];%frequencies for all spectrograms in figure;

% select main folder where we can find subfolder names & excel file
[~,mainpath] = uiputfile('*.*','Select main data folder', 'mainpath.mat');
cd(mainpath)
% Read excel sheet for getting acq #
LoadExcel('allPhotometryData.xlsx')
% batchTableRaw{8,1} = {'2_2'}; % correct for read-in error

for animal_i= [38:49]%[138:141] %85:size(batchTableRaw,1) %[8,25] %
%     try
    clearvars -except batchTableRaw animal_i mainpath params

    % pathlocation
    FolderPath = strcat(mainpath,batchTableRaw{animal_i,2},'/photometry data/',string(batchTableRaw{animal_i,1}),'/',string(batchTableRaw{animal_i,1}),'_',batchTableRaw{animal_i,3});
    % move to folder
    cd(FolderPath)

    % read-in metadata
    % folder = uigetdir('Choose folder with photometry data to be analyzed');
    % cd(folder);
    load('allData2.mat');

    %% Option 1: clean up flourescence (BS' snippet)
    if 1
%     winLen=500;
%     gg1090=1e-5*ones(2, length(rawSig1));
%     for counter=(winLen+1):(length(rawSig1)-winLen)
%       gg1090(:,counter)=prctile(rawSig1(counter+(-winLen:winLen)),[5 95],'all');
%     end
%     ggg=(rawSig1-gg1090(1,:))./(gg1090(2,:)-gg1090(1,:));
%     rr1090=1e-5*ones(2, length(rawSig2));
%     for counter=(winLen+1):(length(rawSig2)-winLen)
%       rr1090(:,counter)=prctile(rawSig2(counter+(-winLen:winLen)),[5 95],'all');
%     end
%     rrr=(rawSig2-rr1090(1,:))./(rr1090(2,:)-rr1090(1,:));
%     figure; plot(ggg,'g');hold on;plot(rrr,'r')
%     save(strcat('detrend1.mat'),'gg1090','ggg','rr1090','rrr');
    

    detrendWindowTime = 60; % in s
    Carriers = lcm(floor(params.samplerate/params.freq1),floor(params.samplerate/params.freq2));
    rawDetrendWindow = Carriers * floor(detrendWindowTime*params.samplerate/Carriers); % in seconds
                                
    green_m=movmean(rawSig1, rawDetrendWindow);
    green_s=movstd(rawSig1, rawDetrendWindow);
    ggg=(rawSig1-green_m)./green_s;
    g_stdZeros=find(green_s==0);
    if ~isempty(g_stdZeros)
        disp('WARNING: Found zeros in standard deviation.  Should not happen.  Setting Infs to 0');
        ggg(g_stdZeros)=0;
    end

    red_m=movmean(rawSig2, rawDetrendWindow);
    red_s=movstd(rawSig2, rawDetrendWindow);
    rrr=(rawSig2-red_m)./red_s;
    r_stdZeros=find(red_s==0);
    if ~isempty(r_stdZeros)
        disp('WARNING: Found zeros in standard deviation.  Should not happen.  Setting Infs to 0');
        rrr(r_stdZeros)=0;
    end
    
    % BS's function to ensure z-score is transformed to normal -- this is
    % for across-day comparison
    ggg_norm=normData(ggg,zCut=0);
    rrr_norm=normData(rrr,zCut=0);
    if ~isempty(tone)
        plotTrialAvgRes(tone,shock,ggg_norm,rrr_norm)
    end
    end

    %% Option 2: Old downsampling and exponential decay fit
    if 0
        params.dsRate= 0.001; % in seconds

        % Determine time points in output data set
        dsTimes = filtTimes(1):params.dsRate:filtTimes(end); 

        % Fit each channel to exponential decay for photobleaching correction
        pbFit1 = fit(filtTimes',filtSig1','exp1'); 
        pbFit2 = fit(filtTimes',filtSig2','exp1');

        % Divide filtered signal by exponential decay to correct for photobleaching
        pbVals1 = (filtSig1' ./ double(pbFit1(filtTimes)))';
        pbVals2 = (filtSig2' ./ double(pbFit2(filtTimes)))';

        % Downsample photobleaching-corrected data
        dsVals1 = interp1(filtTimes,pbVals1,dsTimes,'spline');
        dsVals2 = interp1(filtTimes,pbVals2,dsTimes,'spline');

        % Downsample raw (uncorrected) data
        rawVals1 = interp1(filtTimes,filtSig1',dsTimes,'spline');
        rawVals2 = interp1(filtTimes,filtSig2',dsTimes,'spline');

        % Fs=20;
        % g_0 = running_percentile(dsVals1, floor(Fs*60), 10);
        % g_dff = dsVals1'./g_0 - 1;
        % g_dff = g_dff';

        % figure
        % subplot(3,1,1); hold on
        % plot(rawVals1,'g'); title('low pass filtered spec power downsampled');
        % subplot(3,1,2); hold on
        % plot(rawVals2,'r');
        % subplot(3,1,3); hold on
        % plot(rawVals3,'r');
        % subplot(4,1,4);hold on
        % plot(filtTimes,rawSig3,'g'); title('spec power of excitation light');
        % plot(filtTimes,rawSig4,'r');
    end
    %% Option 3: Sarah's fit? 
    if 0
        startfit=round(find(light==0,1,'first')/20);
        endfit=round(find(light==0,1,'last')/20);

        timeoutloc=[(find(rawSig1(50000:end)<0.8,1,'first')-5000+49000),(find(rawSig1<0.8,1,'last'))+500];
        time=zeros(1,length(startfit:5000:timeoutloc(1)));
        for i=1:length(time)
            [~,b]=min(rawSig1(((i*5000)-5000+startfit):((i*5000)+startfit-1)));
            b=b+((i*5000)-5000+startfit)-1;
            time(i)=b;
        end
        filtSig1fit = fit(time',(rawSig1subtract(time))','exp2');
        h=figure;plot(rawSig1subtract);hold on;plot(filtSig1fit(1:timeoutloc(1)));
        filtSig1fitted=rawSig1subtract;
        filtSig1fitted(1:timeoutloc(1))=filtSig1fitted(1:timeoutloc(1))-filtSig1fit(1:timeoutloc(1))';

        time=zeros(1,length(timeoutloc(2)+10:5000:round((find(light==0,1,'last'))/20)-5000));
        for i=1:length(time)
            [~,b]=min(rawSig1(((i*5000)-5000+10+timeoutloc(2)):((i*5000)+9+timeoutloc(2))));
            b=b+((i*5000)-5000+10+timeoutloc(2))-1;
            time(i)=b;
        end
        filtSig1fit = fit(time',(rawSig1subtract(time))','exp2');
        plot(timeoutloc(2):length(rawSig1),filtSig1fit(timeoutloc(2):length(rawSig1)));
        fitok=input('fit ok???');
        close(h)
        filtSig1fitted(timeoutloc(2):length(rawSig1))=filtSig1fitted(timeoutloc(2):length(rawSig1))-filtSig1fit(timeoutloc(2):length(rawSig1))';
        filtSig1fitted(timeoutloc(1):timeoutloc(2))=1;
        fitfinal_GCaMP=mean(filtSig1fit(round((find(light==0,1,'last'))/20)-5000:round((find(light==0,1,'last'))/20)));
        A=sort(filtSig1fitted(round((find(light==0,1,'first'))/20):round((find(light==0,1,'first'))/20)+6075));
        baseline90_GCaMP=A(round(length(A)*0.9));
        Sig1=filtSig1fitted;
        plot(Sig1);
        fitok=input('Sig1 ok???');
        % close all
    end
    
%% save metadata
% WITHOUT detrend step before demod
%     save('allData.mat','green','red','modgreen','modred','rawSig1','rawSig2','filtSig1','filtSig2','shock','tone','light','filtTimes');
% WITH detrend step before demod
save('allData2_Step2.mat','ggg','rrr','ggg_norm','rrr_norm')
close all
%     end
end