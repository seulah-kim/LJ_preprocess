function [spectTimes, filtSig, sig, spectFreqs, spectAmpVals, returnError] = ...
    spectBS(rawData, params)

% params.winSize2 = 1;
% params.spectSample2 = params.winSize2 ./ 5;

% Reads in:
% 1) rawData -- the raw oscillating photometry signal
% 2) params.rawSampleFreq -- sampling frequency
% 3) freqRange -- 2 element array of frequency ranges to be analyzed
% 4) params -- structure with following elements (and default values used):

% Scott Owen -- 2018-08-12

% Convert spectrogram window size and overlap from time to samples

%disp(['Spectrum window ', num2str(params.spectralWindow ./ params.rawSampleFreq), ' sec; ',...
%    num2str(params.spectralWindow), ' samples at ', num2str(params.rawSampleFreq), ' Hz'])

% Calculate spectrogram

returnError={};
[spectVals,spectFreqs,spectTimes]=spectrogram(rawData,params.spectralWindow,params.spectralWindowOverlap,params.freqRange,params.rawSampleFreq);

% Convert spectrogram to real units
spectAmpVals = double(abs(spectVals));
avgFreqAmps = mean(spectAmpVals,2);

% find frequency with peak power
[~,maxFreqBin]=max(double(avgFreqAmps));
if maxFreqBin~=(length(params.freqRange)+1)/2
    returnError=statusUpdate(returnError, ...
        'WARNING: Spectrogram:  Peak power is not at the carrier frequency');
    returnError=statusUpdate(returnError, ...
        ['   Peak at ' num2str(params.freqRange(maxFreqBin))]);
end

% Calculate signal at each frequency band
sig = mean(abs(spectVals((maxFreqBin-params.inclFreqWin):(maxFreqBin+params.inclFreqWin),:)),1);

filtFreq=1/(spectTimes(2)-spectTimes(1));
if filtFreq<params.lowPassCorner
    filtSig=sig;
else
    % Create low pass filter for final data
    lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',params.lowPassCorner,...
        'PassbandRipple',0.01, 'SampleRate', filtFreq);

    % Low pass filter the signals
    filtSig=filtfilt(lpFilt,double(sig));
end