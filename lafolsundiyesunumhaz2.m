%% Calling ECG Signal
% Open .mat file
clear; clc; close all;
[FileName,PathName] = uigetfile('*.mat','Select mat file');
if FileName==0,
return;
end
MatFile=load(fullfile(PathName,FileName));   %# pass file path as string
Structname = fieldnames(MatFile);
assignin('base', 'ecgsignal', MatFile.(Structname{1}));
fs=1000;
bpm=60;
%% Baseline Drift Removing 
% Remove the moving artifacts
ecg_signal=ecgsignal-mean(ecgsignal);
[c,l]=wavedec(ecg_signal,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
ecgcorrected=ecg_signal-a8;
fl=15; fu=40; % lower and upper cutoff freqs
wcl=2*fl/fs;
wcu=2*fu/fs;
faxes=linspace(-fs/2,fs/2,length(ecgsignal));
plot(faxes,fftshift(abs(fft(ecgsignal))));
axis([0 100 0 4e6])
N=50;
wn=[wcl wcu];
b=fir1(N,wn,'stop',hamming(N+1)); %notch filter with hamming window
iv=zeros(1,N); %initialization vector for all taps to zero
ecgfiltered_last=filter(b,1,ecgcorrected,iv);
ecgfiltered_last = sgolayfilt(ecgfiltered_last,3,81);
% figure; plot(ecgsignal/max(ecgsignal)); title('Raw signal');
% figure; plot(ecgfiltered_last/max(ecgfiltered_last)); title('Filtered signal');
ecgfiltered_last=ecgfiltered_last/max(ecgfiltered_last);

plot(ecgfiltered_last)
hold on
plot(ecgsignal/max(ecgsignal),'r')