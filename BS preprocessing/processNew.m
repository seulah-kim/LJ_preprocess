%%
% this code is for labjack alignment of dual hemisphere, dual color recordings which
% include Lick and LED detection 8 channels

%%
close all
clear all

doPlot=false; % 0 plot nothing, 1 plot basics, 2 plot everything

skipIfExists=false;
onlyRerunExisting=true;

reLoadIfExists=true;
reRunAllIfExists=false;
rePreprocessIfExists=false;
redoNoiseAnalysisIfExists=false;

reAnalyzeBehavaiorIfExists_table=false;
reAnalyzeBehavaiorIfExists_matrices=true;

reAnalyzeCCsIfexists=true;

if reRunAllIfExists
    rePreprocessIfExists=true;
    redoNoiseAnalysisIfExists=true;
    reAnalyzeBehavaiorIfExists_table=true;
    reAnalyzeBehavaiorIfExists_matrices=true;
    reAnalyzeCCsIfexists=true;
end

% put here the parameters that you want to change on reAnalysis
overwriteParams=struct;
overwriteParams.ptsKeep_after=100;

searchFolders=1;

% saveFolder='/Users/bernardosabatini/Dropbox (HMS)/2ABT_data_bernardo/test'
% photomFolder='/Users/bernardosabatini/Dropbox (HMS)/2ABT_data_bernardo/photometry'
% dateList={'11222021'}
% mouseList={'WT63'}

saveFolder='/Volumes/BS Office/Dropbox (HMS)/2ABT_data_bernardo/new_analysis/';
photomFolder='/Volumes/Neurobio/MICROSCOPE/Lynne/2ABT/1_photom';

%% if searchFolders is set, then move through the directory structure and analyze everything
%   start at the levels of the date folders
if searchFolders
    cd(photomFolder)
    D=flip(dir()); % this is where you can select dates.  Enter dir('01282021*'), for example
    
    dCounter=1;
    dateList={};
    for dfCounter=1:length(D)
        if length(D(dfCounter).name)==8
            dateList{dCounter}=D(dfCounter).name;
            dCounter=dCounter+1;
        end
    end
end

%% Mice we know we don't want to analyze
skipMiceList={'S1291', 'S1292', 'S1293', 'S1294'};


%% Loop through all the dates
for dateC=dateList
    date=dateC{1};
    %    dateFolder=sprintf('/Volumes/BS Office/Dropbox (HMS)/2ABT_data_bernardo/simultaneous recordings/photometry recording/%s/', date);
    
    %% If searchFolders is set, now gather all the mouse by mouse folders
    if searchFolders
        dateFolder=sprintf('/Volumes/Neurobio/MICROSCOPE/Lynne/2ABT/1_photom/%s/', date);
        cd(dateFolder)
        disp(['Scanning ' dateFolder])
        D=dir(); %[dir('WT54*')' dir('WT55*')' dir('WT56*')']';     % this is where you can select mice.  Enter dir('WT63*'), for example
        dCounter=1;
        mouseList={};
        for dfCounter=1:length(D)
            if length(D(dfCounter).name)>2 && ~any(strcmp(D(dfCounter).name, skipMiceList))
                mouseList{dCounter}=D(dfCounter).name;
                dCounter=dCounter+1;
            end
        end
    end
    
    %% Loop through all the mice
    for mouseC=mouseList
        mouse=mouseC{1};
        
        params=struct;
        params.notes={};
        cd(saveFolder)
        
        disp(' ')
        
        if skipIfExists && exist(['processed_' mouse '_' date '.mat'], 'file')
            disp([mouse ' ' date ' already analyzed. Skipping. Change flag skipIfExists to reanalyze'])
        else
            %        folder1 = sprintf('/Volumes/BS Office/Dropbox (HMS)/2ABT_data_bernardo/simultaneous recordings/photometry recording/%s/%s', mouse, date); %photometry
            %        folder2 = sprintf('/Volumes/BS Office/Dropbox (HMS)/2ABT_data_bernardo/simultaneous recordings/behavior recording/%s/%s', mouse, date); %behavior
            %        folder1=sprintf('/Users/bernardosabatini/Dropbox (HMS)/2ABT_data_bernardo/photometry/%s/%s', date, mouse)
            %        folder2=folder1;
            
            disp(['Working on ' mouse ' ' date ]);
            
            if reLoadIfExists && exist(['processed_' mouse '_' date '.mat'], 'file')
                disp('RELOADING')
                temp= load(['processed_' mouse '_' date '.mat']);
                processed=temp.(['processed_' mouse '_' date]);
                clear temp
                params=processed.params;
                folder1=params.photometry_folder;
                folder2=params.behavior_folder;
                reloaded=1;
                overwriteParameters
            else
                folder1 = sprintf('/Volumes/Neurobio/MICROSCOPE/Lynne/2ABT/1_photom/%s/%s', date, mouse); %photometry
                folder2 = sprintf('/Volumes/Neurobio/MICROSCOPE/Lynne/2ABT/1_behavior/%s/%s', date, mouse); %behavior
                reloaded=0;
            end
            
            if ~exist(folder1, 'dir') || ~exist(folder2, 'dir')
                disp('   A data folder is missing. Skipping...');
            else
                %% check if the behavior data exists before progressing
                pokeFile='';
                statsFile='';
                cd(folder2);
                files = dir('.');
                for i=1:length(files)
                    if contains(files(i).name, 'pokeHistory')
                        pokeFile = files(i).name;
                    elseif regexp(files(i).name, 'stats.*\.mat')
                        statsFile = files(i).name;
                    end
                end
                
                if exist(pokeFile, 'file') && exist(statsFile, 'file')
                    %% some basic parameters
                    if (~reloaded && ~onlyRerunExisting) || (reloaded && rePreprocessIfExists)             
                        if ~reloaded
                            params.mouse=mouse;
                            params.date=date;
                            params.photometry_folder=folder1;
                            params.behavior_folder=folder2;
                            
                            params.rawSampleFreq=2000;
                            params.numChannels=14; % number of channels
                        else
                            overwriteParameters;
                        end
                        
                        %% Concatenate data files into an array
                        %data array has the info for each channel that was collected
                        cd(folder1);
                        
                        D=dir('Raw_*.mat');
                        if ~isempty(D)
                            filename={D.name};
                            load(filename{1});
                            nPtsPerTemp=length(temp);
                            disp(['   Loaded first section.  # data points: ' num2str(nPtsPerTemp)]);
                            
                            is13ChanData=0;
                            if nPtsPerTemp==26000
                                disp('   Seems like 13 channel data')
                                is13ChanData=1;
                            elseif nPtsPerTemp==28000
                                disp('   Seems like 14 channel data')
                                is13ChanData=0;
                            else
                                disp('   Channel structure is unknown')
                            end
                            
                            output=zeros(1,(length(filename)*(nPtsPerTemp)));
                            
                            nChansAssume=floor(nPtsPerTemp/params.rawSampleFreq);
                            if nChansAssume==nPtsPerTemp/params.rawSampleFreq
                                params.notes=statusUpdate(params.notes, ...
                                    ['Structuring as ' num2str(nChansAssume) ' channels']);
                                for i=1:length(D)
                                    if i>1
                                        load(filename{i});
                                    end
                                    output(((i-1)*(nPtsPerTemp)+1):(i*(nPtsPerTemp)))=temp;
                                end
                                if (length(output)-floor(length(output)/nChansAssume)*nChansAssume)>eps
                                    params.notes=statusUpdate(params.notes, ...
                                        'WARNING: output is not integer multiple of channel length. TRIMMING...');
                                end
                                clear temp
                            else
                                output=[];
                                params.notes=statusUpdate(params.notes, ...
                                    'WARNING: Unclear how many channels are present. Skipping...');
                            end
                        else
                            output=[];
                            params.notes=statusUpdate(params.notes, ...
                                'WARNING: Directory has no RAW files. Skipping...');
                        end
                        
                        if length(output)>10000 % then continue.  Otherwise, what is this?
                            % reshape loaded data
                            params.samplesPerChannel=floor(length(output)/nChansAssume);
                            output=reshape(output(1:(params.samplesPerChannel*nChansAssume)), nChansAssume, params.samplesPerChannel);
                            if nChansAssume<params.numChannels
                                output(end+1:params.numChannels, :)=0;
                            end
                            
                            %% plot all raw data
                            
                            if is13ChanData
                                params.channelDefs={...
                                    'AIN0', ...                 % Ch 1
                                    'AIN1', ...                 % Ch 2
                                    'AIN2', ...                 % Ch 3
                                    'AIN3', ...                 % Ch 4
                                    'AIN4', ...                 % Ch 5
                                    'AIN5', ...                 % Ch 6
                                    'Centerport', ...           % Ch 7
                                    'Rightport', ...            % Ch 8
                                    'Leftport', ...             % Ch 9
                                    'Left Lick', ...            % Ch 10
                                    'Right Lick', ...           % Ch 11
                                    'Center LED', ...           % Ch 12
                                    'Left LED', ...             % Ch 13
                                    'Zeros' ...                 % Ch 14
                                    };
                            else
                                params.channelDefs={...
                                    'AIN0', ...                 % Ch 1
                                    'AIN1', ...                 % Ch 2
                                    'AIN2', ...                 % Ch 3
                                    'AIN3', ...                 % Ch 4
                                    'AIN4', ...                 % Ch 5
                                    'AIN5', ...                 % Ch 6
                                    'Centerport', ...           % Ch 7
                                    'Rightport', ...            % Ch 8
                                    'Leftport', ...             % Ch 9
                                    'Left Lick', ...            % Ch 10
                                    'Right Lick', ...           % Ch 11
                                    'Center LED', ...           % Ch 12
                                    'Left LED' ...              % Ch 13
                                    'Laser' ...                 % Ch 14
                                    };
                            end

                            centerPortChannel=find(strcmp(params.channelDefs, 'Centerport'));

                            if doPlot
                                figure;
                                for channel=1:params.numChannels
                                    subplot(params.numChannels, 1, channel);
                                    plot(output(channel, :));
                                    title(params.channelDefs{channel});
                                end
                            end
                            
                            %% Set up info on the channels
                            
                            params.channelNames={...
                                'green r', ...      % channel 1
                                'red r', ...        % channel 2
                                'gCarrier', ...     % channel 3
                                'rCarrier', ...     % channel 4
                                'green l', ...      % channel 5
                                'red l' ...         % channel 6
                                };
                            params.channelNames{15}='corrected red r'; % channel 15
                            params.channelNames{16}='corrected red l'; % channel 16
                            
                            %% Set up spectrogram parameters
                            % Instructions
                            % 1) Set min and max values for params.freqRange1 and params.freqRange2 so
                            %    each range spans the target oscillation frequency of that channel, and
                            %    (ideally) excludes the oscillation on the opposite channel
                            % 2) Start the script
                            % 3) Most useful output variables are:
                            %    dsTimes             -- the time stamps of the output data
                            %    rawVals1 & rawVals2 -- output data before photobleaching correction
                            %    dsVals1 & dsVals2   -- photobleaching corrected data for channels 1 and 2
                            %    syncInTimes         -- timestamps of TTL events on SyncIn channel
                            %    syncOutTimes        -- timestamps of TTL events on SyncOut channel
                            
                            params.finalSampleFreq=params.rawSampleFreq/(12*9);
                            params.finalTimeStep=1/params.finalSampleFreq;
                            params.inclFreqWin = 4; % Number of frequency bins to average (on either side of peak freq)
                            params.freqStep = 1; % Step sfize in Hz for freqency calculations
                            params.freqStepWidth = 7;
                            params.detrendWindowTime = 60; % in seconds
                            params.dropFirstDetrendWindow=1; % drop the signals before the end of the first detrend time window
                            params.lowPassCorner=100;
                            params.ptsKeep_before=40;
                            params.ptsKeep_after=60;
                            params.spectralWindow=2*9*12;
                            params.spectralWindowOverlap=params.spectralWindow/2;
                            params.carrierChannel=[3 4 3 4 3 4]; % which channel has the carrier for each fluorescence channel
                            params.carrierChannel(15:16)=4;
                            
                            if reloaded
                                overwriteParameters;
                            end

                            %% Clean up extra zeros, if there are any
                            outputChannelMean=mean(output,1);
                            outputMeanZeros=find(outputChannelMean==0);
                            if ~isempty(outputMeanZeros)
                                lastKeep=min(outputMeanZeros);
                                params.notes=statusUpdate(params.notes, ...
                                    ['WARNING: Found all channel zeros in raw data.  Truncating to ' ...
                                    num2str(lastKeep) ' points']);
                                output=output(:, 1:lastKeep);
                            end

                            clear outputChannelMean outputMeanZeros
                            
                            %% correct cross talk
                            params.mdl=cell(1,6);
                            mdl = fitlm(output(1,:),output(2,:));
                            output(15,:)=output(2,:)-mdl.Coefficients.Estimate(2)*output(1,:);
                            params.mdl{2}=mdl;
                            
                            mdl = fitlm(output(5,:),output(6,:));
                            output(16,:)=output(6,:)-mdl.Coefficients.Estimate(2)*output(5,:);
                            params.mdl{6}=mdl;
                            
                            params.numChannels=params.numChannels+2;
                            
                            %% Set up
                            channelsToAnalze=[1 2 5 6];% 15 16];
                            channelsToGetCarrierFreq=[1 2 3 4 5 6];% 15 16];
                            
                            %% Extract carrier frequencies
                            params.measuredCarrierFreq=zeros(1,6);
                            params.carrierPtsPerCycle=zeros(1,6);
                            params.setCarrierFreq=[167 223 167 223 167 223];
                            
                            pointsToProcess=2^14;
                            for channel=channelsToGetCarrierFreq
                                ChTitleString=[' Ch' num2str(channel) ' (' params.channelNames{channel} ') '];
                                Y=fft(normalize(output(channel, 1:pointsToProcess)));
                                P2 = abs(Y/pointsToProcess);
                                P1 = P2(1:pointsToProcess/2+1);
                                P1(2:end-1) = 2*P1(2:end-1);
                                f = params.rawSampleFreq*(0:(pointsToProcess/2))/pointsToProcess;
                                if doPlot==2
                                    figure; plot(f,P1)
                                    title(['FFT of Ch' num2str(channel) ': ' params.channelNames{channel} ]);
                                end
                                [~, maxFindex]=max(P1);
                                candFreq=f(maxFindex);
                                params.carrierPtsPerCycle(channel)=floor(params.rawSampleFreq/candFreq);
                                params.measuredCarrierFreq(channel)=params.rawSampleFreq/params.carrierPtsPerCycle(channel);
                                params.notes=statusUpdate(params.notes, ...
                                    ['For Ch' num2str(channel) ' (' params.channelNames{channel} ...
                                    ') the dominant frequency is ' num2str(params.measuredCarrierFreq(channel))]);
                            end
                            
                            clear Y P1 P2 f
                            params.lcmCarriers=lcm(params.carrierPtsPerCycle(3), params.carrierPtsPerCycle(4));
                            params.notes=statusUpdate(params.notes, ...
                                ['   Assuming carriers are in channels 3 and 4, the recommended minimum analysis unit is ' ...
                                num2str(params.lcmCarriers) ' points']);
                            
                            
                            %% Process the data
                            
                            params.useMeasuredFreq=1;
                            params.finalSamples=0;
                            params.totalDownSample=params.rawSampleFreq*params.finalTimeStep;

                            if reloaded
                                overwriteParameters;
                            end
                            
                            params.notes=statusUpdate(params.notes, ...
                                'Processing flourescence signals');
                            firstRun=1;
                            %%
                            for channel=channelsToAnalze
                                chString=['Ch' num2str(channel) ' (' params.channelNames{channel} ')'];
                                %                        disp(['*** processing ' chString]);
                                
                                if (params.carrierPtsPerCycle(channel) - ...
                                        params.carrierPtsPerCycle(params.carrierChannel(channel))) > eps
                                    params.notes=statusUpdate(params.notes, ...
                                        ['WARNING: ' chString ' carrier frequency is wrong. Using designated carrier']);
                                    params.notes=statusUpdate(params.notes, ...
                                        '   Assuming large cross talk.  Try examining Channels 15 and 16 for red..');
                                    freq=params.measuredCarrierFreq(params.carrierChannel(channel));
                                else
                                    if params.useMeasuredFreq
                                        freq=params.measuredCarrierFreq(channel);
                                    else
                                        freq=params.setCarrierFreq(channel);
                                    end
                                end
                                
                                params.freqRange = freq + ...
                                    params.freqStep*...
                                    (-params.freqStepWidth:params.freqStepWidth); % define a frequency band around peak
                                
                                flRaw=output(channel, :);
                                
                                % here we force the detrending window to be
                                % a multiple of the LCM of the carrier
                                % signals periods
                                params.rawDetrendWindow = ...
                                    params.lcmCarriers ...
                                    * floor(params.detrendWindowTime*params.rawSampleFreq/params.lcmCarriers); % in seconds
                                
                                flRaw_m=movmean(flRaw, params.rawDetrendWindow);
                                flRaw_s=movstd(flRaw, params.rawDetrendWindow);
                                flSignal=(flRaw-flRaw_m)./flRaw_s;
                                stdZeros=find(flRaw_s==0);

                                % keep some copy of the mean from
                                % detrending to calculate dF/F
                                blockLength=params.totalDownSample;
                                nBlocks=floor(length(flRaw_m)/blockLength);
                                flDownSample_f0=mean(...
                                    reshape(flRaw_m(1:(blockLength*nBlocks)), blockLength, nBlocks),...
                                    1);

                                if ~isempty(stdZeros)
                                    params.notes=statusUpdate(params.notes, ...
                                        'WARNING: Found zeros in standard deviation.  Should not happen.  Setting Infs to 0');
                                    flSignal(stdZeros)=0;
                                end
                                
                                % Calculate spectrogram channel 1
                                try
                                    [filtTimes, filtSignal, rawSignal, spectFreqs, spectAmpVals, returnError] = spectBS(...
                                        flSignal, ...
                                        params...
                                        );
                                    params.notes=statusUpdate(params.notes, returnError);
                                catch ME
                                    params.notes=statusUpdate(params.notes, ...
                                        'WARNING: Error in spectrogram!');
                                    filtSignal=[];
                                end
                                
                                params.freqRange=[];
                                
                                if ~isempty(filtSignal)
                                    filtdt=filtTimes(2)-filtTimes(1);
                                    filtt0=filtTimes(1);
                                    
                                    blockLength=params.finalTimeStep/filtdt;
                                    if (floor(blockLength)-blockLength)>eps
                                        error('     Spectral filter sample frequency is not integer of desired');
                                    else
                                        blockLength=floor(blockLength);
                                    end
                                    
                                    nBlocks=floor(length(filtSignal)/blockLength);
                                    dSignal=mean(...
                                        reshape(filtSignal(1:(blockLength*nBlocks)), blockLength, nBlocks),...
                                        1);
                                    
                                    params.signalDetrendWindow = ...
                                        floor(params.detrendWindowTime*params.finalSampleFreq);
                                    
                                    pSignal=rollingZ(dSignal, params.signalDetrendWindow);
                                    
                                    params.finalSamples=length(pSignal);
                                    t0=filtTimes(1);
                                    
                                    params.final_t0(channel)=t0;
                                    params.final_samples(channel)=params.finalSamples;
                                    
                                    % Determine time points in output data set
                                    flSampleTimes = params.final_t0(channel)+...
                                        params.finalTimeStep*(1:params.finalSamples);
                                    
                                    if firstRun
                                        processed=struct;
                                        processed.signals=cell(1, params.numChannels);
%                                         processed.rawf0=cell(1, params.numChannels);
%                                         processed.finalf0=cell(1, params.numChannels);

                                        firstRun=0;
                                    end
                                    
                                    % the detrended signal
                                    processed.signals{channel}=pSignal;

%                                     processed.rawf0{channel}=flDownSample_f0(1:length(pSignal));
%                                     processed.finalf0{channel}=movmean(dSignal, params.signalDetrendWindow);
                                end
                            end
                        end
                    end
                    
                    if ~isfield(params, 'finalSamples')
                        params.finalSamples=0;
                    end
                    
                    if params.finalSamples>0
                        if ~reloaded || (reloaded && redoNoiseAnalysisIfExists)
                            for ccc=1:length(processed.signals)
                                if ~isempty(processed.signals{ccc})
                                    processed.signalMoments(ccc, 1)=mean(processed.signals{ccc});
                                    processed.signalMoments(ccc, 2)=var(processed.signals{ccc});
                                    processed.signalMoments(ccc, 3)=skewness(processed.signals{ccc}, 1);
                                    processed.signalMoments(ccc, 4)=kurtosis(processed.signals{ccc}, 1);
                                    processed.signalMoments(ccc, 5)=skewness(processed.signals{ccc}, 0); % correct bias
                                    processed.signalMoments(ccc, 6)=kurtosis(processed.signals{ccc}, 0); % correct bias
                                end
                            end
                        end

                        if ~reloaded || (reloaded && rePreprocessIfExists)
                            params.behaviorRange=7:14;
                            params.behaviorRangeN=length(params.behaviorRange);
                            params.notes=statusUpdate(params.notes, ...
                                'Down sampling behavior event data');
                            params.totalDownSample=params.rawSampleFreq*params.finalTimeStep;
                            outputSignalRange=(1:params.finalSamples*params.totalDownSample); %+t0*params.rawSampleFreq
                            
                            processed.behavior.downSampled=...
                                squeeze(sum(...
                                reshape(...
                                output(params.behaviorRange,outputSignalRange), ...
                                params.behaviorRangeN, params.totalDownSample, params.finalSamples), ...
                                2));
                            
                            processed.behavior.risingEdge=...
                                squeeze(sum(...
                                reshape(...
                                [diff(output(params.behaviorRange,outputSignalRange), 1, 2) zeros(params.behaviorRangeN, 1)]==1, ...
                                params.behaviorRangeN, params.totalDownSample, params.finalSamples), ...
                                2))>0;
                            
                            processed.behavior.fallingEdge=...
                                squeeze(sum(...
                                reshape(...
                                [diff(output(params.behaviorRange,outputSignalRange), 1, 2) zeros(params.behaviorRangeN, 1)]==-1, ...
                                params.behaviorRangeN, params.totalDownSample, params.finalSamples), ...
                                2))>0;
                            
                            processed.behavior.occupance=...
                                squeeze(sum(...
                                reshape(...
                                output(params.behaviorRange,outputSignalRange), ...
                                params.behaviorRangeN, params.totalDownSample, params.finalSamples), ...
                                2))>0;
                        end
                        
                        behaviorErrorFlag=false;
                        if ~reloaded || (reloaded && (reAnalyzeBehavaiorIfExists_table || reAnalyzeBehavaiorIfExists_matrices))
                            if reloaded
                                overwriteParameters;
                            end
                            processBehavior;
                        end
                        
                        if ~behaviorErrorFlag
                            fullSaveFile=fullfile(saveFolder, ['processed_' mouse '_' date '.mat']);
                            useCurrent=1;
                            processed.params=params;
                            if ~reloaded || (reloaded && reAnalyzeCCsIfexists)
                                anaCrossCor
                            end
                            params.notes=statusUpdate(params.notes, ...
                                ['Saving to ' fullSaveFile]);
                            assignin('base', ['processed_'  mouse '_' date], processed);
                            save(fullSaveFile, ['processed_'  mouse '_' date])
                            eval(['clear processed_'  mouse '_' date]);
                            clear processed
                            clear params
                            clear output
                        else
                            params.notes=statusUpdate(params.notes, ...
                                'WARNING: Analysis failed');
                        end
                    else
                        params.notes=statusUpdate(params.notes, ...
                            'WARNING: Too little data. Skipping');
                    end
                else
                    params.notes=statusUpdate(params.notes, ...
                        'WARNING: A behavior file is missing. Skipping...');
                end
            end
        end
    end
end

