% PreprocessLJStep1.m
% This file reads Raw_****.mat generated by LabJack. Datapoints will be sorted based on scanning order into photometry and behavior data outputs
% The photometry signals are demodulated using spectrogram analysis.
% Spectrogram plot is generated

close all
% set parameters
params.samplerate=2052; % Hz
params.spectWindow=200; % window size for frequency calculation in spectrogram (number of samples) =97.5ms
params.spectOverlap=180; % overlap between windows in spectrogram; new calculation every 20 samples.data points thus downsampled to 102.6Hz.
% params.freqRange1 = [163.5:5:178.5]; % frequencies for channel 1 spectrogram (target 171 Hz) 
params.freqRange1 = [166,171,176]; % frequencies for channel 1 spectrogram (target 171 Hz) 
params.freq1=171;
params.freqRange2 = [223,228,233]; % frequencies for channel 2 spectrogram (target 228 Hz)
% params.freqRange2 = [220.5:5:235.5]; % frequencies for channel 2 spectrogram (target 228 Hz)
params.freq2=228;
params.useFreqRange=[0:5:500];%frequencies for all spectrograms in figure;

% folder = uigetdir('','Choose folder to save photometry data');
% cd(folder);

% select main folder where we can find subfolder names & excel file
[~,mainpath] = uiputfile('*.*','Select main data folder', 'mainpath.mat');
cd(mainpath)
% Read excel sheet for getting acq #
LoadExcel('allPhotometryData.xlsx')
% batchTableRaw{8,1} = {'2_2'}; % correct for read-in error

for animal_i= [67:73,91:97,113:119,135:141] %1:size(batchTableRaw,1)
     
    clearvars -except batchTableRaw animal_i mainpath params

    % pathlocation
    FolderPath = strcat(mainpath,batchTableRaw{animal_i,2},'/photometry data/',string(batchTableRaw{animal_i,1}),'/',string(batchTableRaw{animal_i,1}),'_',batchTableRaw{animal_i,3});
    % move to folder
    cd(FolderPath)

    D=dir('Raw_*.mat');
    filename={D.name};
    load(filename{1});
    numChannels=length(temp)/params.samplerate;
    newfileIdx = length(D);
    dataArray=zeros(1,(newfileIdx*length(temp)));
    for i=1:newfileIdx
        load(filename{i});
        dataArray(((i-1)*length(temp)+1):(i*length(temp)))=temp;
    end

    output=dataArray;
    totalLen = length(dataArray);

    % numChannels=6; % number of channels (x)
    t_s = (1/params.samplerate)*(1:round(totalLen/numChannels)); %time base in seconds for x axis of figure

    % In the order of scanning
    Ch1=find(mod(1:totalLen,numChannels)==1);   % photodetector #1 -green
    Ch2=find(mod(1:totalLen,numChannels)==2);   % photodetector #2 -red
    Ch3=find(mod(1:totalLen,numChannels)==3);  %copy of 488 modulation
    Ch4=find(mod(1:totalLen,numChannels)==4);   %copy of 560 modulation
    Ch5=find(mod(1:totalLen,numChannels)==5);   %led
    Ch6=find(mod(1:totalLen,numChannels)==6);   %tone
    Ch7=find(mod(1:totalLen,numChannels)==0);   %shock

    % Make green and red arrays and demodulate photodiode signals  
    green= output(Ch1);
    red= output(Ch2);
    modgreen=output(Ch3);
    modred=output(Ch4);
    light=output(Ch5);
    tone=output(Ch6);
    shock=output(Ch7);
    if strcmp(batchTableRaw{animal_i,3},'OF')
        light=output(Ch7);
    end

    % Remove artifacts -- 1/12/21 SAK added
    idx=find(light>1);
    for i=idx
        green(i)=mean([green(i-1),green(i+1)]);
        red(i)=mean([red(i-1),red(i+1)]);
        modgreen(i)=mean([modgreen(i-1),modgreen(i+1)]);
        modred(i)=mean([modred(i-1),modred(i+1)]);
        light(i)=mean([light(i-1),light(i+1)]);
        if ~strcmp(batchTableRaw{animal_i,3},'OF')
            tone(i)=mean([tone(i-1),tone(i+1)]);
            shock(i)=mean([shock(i-1),shock(i+1)]);
        end
    end

    figure('Position',[440 126 571 672]);hold on
    subplot(7,1,1); plot(green);title('green raw');axis tight
    subplot(7,1,2); plot(red);title('red raw');axis tight
    subplot(7,1,3); plot(modgreen);title('green modulation');axis tight
    subplot(7,1,4); plot(modred);title('red modulation');axis tight
    subplot(7,1,5); plot(shock);title('shock');axis tight
    subplot(7,1,6); plot(tone);title('tone');axis tight
    subplot(7,1,7); plot(light);title('light');axis tight
    saveas(gcf,'AllChannelsRaw')
    %% SAK added 8.26.21
    % crop_idx=1000;
    % green= green(crop_idx:end);
    % red= red(crop_idx:end);
    % modgreen=modgreen(crop_idx:end);
    % modred=modred(crop_idx:end);
    % 
    % t_s = t_s(crop_idx:end);
    % %plot ttls
    % light = light(crop_idx:end);
    % tone = tone(crop_idx:end);
    % shock = shock(crop_idx:end);
    %% detrend before demodulation
%     winLen=500;
%     winLen=50000;
%     gg1090=1e-5*ones(2, length(green));
%     for counter=(winLen+1):(length(green)-winLen)
%       gg1090(:,counter)=prctile(green(counter+(-winLen:winLen)),[5 95],'all');
%     end
%     ggg=(green-gg1090(1,:))./(gg1090(2,:)-gg1090(1,:));
%     rr1090=1e-5*ones(2, length(red));
%     for counter=(winLen+1):(length(red)-winLen)
%       rr1090(:,counter)=prctile(red(counter+(-winLen:winLen)),[5 95],'all');
%     end
%     rrr=(red-rr1090(1,:))./(rr1090(2,:)-rr1090(1,:));

if 1   % moving detrend
    detrendWindowTime = 60; % in s
    Carriers = lcm(floor(params.samplerate/params.freq1),floor(params.samplerate/params.freq2));
    rawDetrendWindow = Carriers * floor(detrendWindowTime*params.samplerate/Carriers); % in seconds
                                
    green_m=movmean(green, rawDetrendWindow);
    green_s=movstd(green, rawDetrendWindow);
    ggg=(green-green_m)./green_s;
    g_stdZeros=find(green_s==0);
    if ~isempty(g_stdZeros)
        disp('WARNING: Found zeros in standard deviation.  Should not happen.  Setting Infs to 0');
        ggg(g_stdZeros)=0;
    end

    red_m=movmean(red, rawDetrendWindow);
    red_s=movstd(red, rawDetrendWindow);
    rrr=(red-red_m)./red_s;
    r_stdZeros=find(red_s==0);
    if ~isempty(r_stdZeros)
        disp('WARNING: Found zeros in standard deviation.  Should not happen.  Setting Infs to 0');
        rrr(r_stdZeros)=0;
    end
    
    figure; plot(ggg,'g');hold on;plot(rrr,'r')
    save(strcat('detrend_raw.mat'),'ggg','green_m','green_s','rrr','red_m','red_s');
end
    %%
    % params.spectSample = 0.01; % Step size for spectrogram (sec)
    params.filtCut = 100/(params.spectWindow-params.spectOverlap); % Cut off frequency of 5Hz for low pass filter of processed data
    params.dsRate = 0.05; % Time steps for down-sampling (seconds) (average of every 100 samples)
    
    % SAK added 3.11.22 to try detrend before demod
    green = ggg;
    red = rrr;
    
    % Convert spectrogram window size and overlap from time to samples
    % spectWindow = 2.^nextpow2(samplerate .* params.winSize);
    % spectOverlap = ceil(spectWindow - (spectWindow .* (params.spectSample ./ params.winSize)));
    disp(['Spectrum window ', num2str(params.spectWindow ./ params.samplerate), ' sec; ',...
        num2str(params.spectWindow), ' samples at ', num2str(params.samplerate), ' Hz with ',num2str(params.spectOverlap),' samples overlap'])
    % Create low pass filter for final data
    lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',params.filtCut,...
        'PassbandRipple',0.01, 'Samplerate',params.samplerate/(params.spectWindow-params.spectOverlap));

    % Calculate spectrogram channel 1
    [spectVals1,spectFreqs1,filtTimes]=spectrogram(green,params.spectWindow,params.spectOverlap,params.freqRange1,params.samplerate);
    figure;subplot(2,3,1);hold on; title('green');spectrogram(green,params.spectWindow,params.spectOverlap,params.useFreqRange,params.samplerate);
    rawSig1 = mean(abs(spectVals1),1);
    filtSig1 = filtfilt(lpFilt,double(rawSig1));% Low pass filter the signals
    1

    % Calculate spectrogram channel 2 modulated by green LED
    [spectVals2,spectFreqs2,filtTimes]=spectrogram(red,params.spectWindow,params.spectOverlap,params.freqRange2,params.samplerate);
    subplot(2,3,2);hold on; title('red');spectrogram(red,params.spectWindow,params.spectOverlap,params.useFreqRange,params.samplerate);
    rawSig2 = mean(abs(spectVals2),1);
    filtSig2 = filtfilt(lpFilt,double(rawSig2));% Low pass filter the signals
    2
    % Calculate spectrogram channel 2 modulated by blue LED
    [spectVals3,spectFreqs3,filtTimes]=spectrogram(red,params.spectWindow,params.spectOverlap,params.freqRange1,params.samplerate);
    subplot(2,3,3);hold on; title('red');spectrogram(red,params.spectWindow,params.spectOverlap,params.useFreqRange,params.samplerate);
    rawSig3 = mean(abs(spectVals3),1);
    filtSig3 = filtfilt(lpFilt,double(rawSig3));% Low pass filter the signals
    3
    % Calculate spectrogram channel 3
    [spectVals4,spectFreqs4,filtTimes]=spectrogram(modgreen,params.spectWindow,params.spectOverlap,params.freqRange1,params.samplerate);
    subplot(2,2,3);hold on; title('modgreen');spectrogram(modgreen,params.spectWindow,params.spectOverlap,params.useFreqRange,params.samplerate);
    rawSig4 = mean(abs(spectVals4),1);
    4

    % Calculate spectrogram channel 4
    [spectVals5,spectFreqs5,filtTimes]=spectrogram(modred,params.spectWindow,params.spectOverlap,params.freqRange2,params.samplerate);
    subplot(2,2,4);hold on; title('modred');spectrogram(modred,params.spectWindow,params.spectOverlap,params.useFreqRange,params.samplerate);
    rawSig5 = mean(abs(spectVals5),1);
    5
    saveas(gcf,'spectrogram')

    figure;set(gcf,'color','w')
    subplot(3,1,1); hold on
    plot(filtTimes,rawSig1,'g'); title('green demod with blue freq');
    subplot(3,1,2); hold on
    plot(filtTimes,rawSig2,'r'); title('red demod with amber freq');
    subplot(3,1,3); hold on
    plot(filtTimes,rawSig3,'r'); title('red demod with blue freq');
    saveas(gcf,'demodsubplots')
    figure;set(gcf,'color','w')
    subplot(3,1,1); hold on
    plot(filtTimes,filtSig1,'g'); title('LP filtered:green demod with blue freq');
    subplot(3,1,2); hold on
    plot(filtTimes,filtSig2,'r'); title('LP filtered:red demod with amber freq');
    subplot(3,1,3); hold on
    plot(filtTimes,filtSig3,'r'); title('LP filtered:red demod with blue freq');
    saveas(gcf,'demodsubplots_LPfiltered')

    %% save metadata
    % WITHOUT detrend step before demod
%     save('allData.mat','green','red','modgreen','modred','rawSig1','rawSig2','filtSig1','filtSig2','shock','tone','light','filtTimes');
    % WITH detrend step before demod
    save('allData2.mat','green','red','modgreen','modred','rawSig1','rawSig2','filtSig1','filtSig2','shock','tone','light','filtTimes');
    close all
%     run Analysis_sound_shock_photometry_Seulah_FC_SAK082621
end