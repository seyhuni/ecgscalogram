close all; clear all; clc;
% Saðlýklý
clear; clc; close all;
[FileName,PathName] = uigetfile('*.mat','Select mat file');
if FileName==0,
return;
end
MatFile=load(fullfile(PathName,FileName));   %# pass file path as string
Structname = fieldnames(MatFile);
assignin('base', 'ecgsignal', MatFile.(Structname{1}));
fs=1000;
%% Baseline Drift Removing 
% Remove the moving artifacts
ecg_signal=ecgsignal-mean(ecgsignal);
[c,l]=wavedec(ecg_signal,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
ecgcorrected=ecg_signal-a8;
%% Adaptive Filtering
fl=40; fu=60; % lower and upper cutoff freqs
wcl=2*fl/fs;
wcu=2*fu/fs;
faxes=linspace(-fs/2,fs/2,length(ecgsignal));
% plot(faxes,fftshift(abs(fft(ecgsignal))));
% axis([0 100 0 4e6])
N=100;
wn=[wcl wcu];
b=fir1(N,wn,'stop',hann(N+1)); %notch filter with hamming window
[h,w]= freqz(b,1,256);
j=1:length(ecgcorrected);
iv=zeros(1,N); %initialization vector for all taps to zero
ecgfiltered_last=filter(b,1,ecgcorrected,iv);
ecgfiltered_last=ecgfiltered_last/max(ecgfiltered_last);
plot(ecgfiltered_last(1:3000))
b=ecgfiltered_last(1086:1271);
save a b;
clear;
wt=cwt(b,1:256,'db6','scal');
figure
sc1=wscalogram('image',wt);
close all;
plot(b);
p=sum(sum(sc1(:,1:35)))
qrs=sum(sum(sc1(:,50:91)))
t=sum(sum(sc1(:,92:130)))
close all;
