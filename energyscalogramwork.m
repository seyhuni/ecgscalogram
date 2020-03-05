tic
%% Calling ECG Signal
% Open .mat file
clear;clc;close all;
[FileName,PathName] = uigetfile('*.mat','Select mat file');
if FileName==0,
return;
end
MatFile=load(fullfile(PathName,FileName));   %# pass file path as string
Structname = fieldnames(MatFile);
assignin('base', 'ecg1', MatFile.(Structname{1}));
%% Calling ECG Signal
% Open .mat file
[FileName,PathName] = uigetfile('*.mat','Select mat file');
if FileName==0,
return;
end
MatFile=load(fullfile(PathName,FileName));   %# pass file path as string
Structname = fieldnames(MatFile);
assignin('base', 'ecg2', MatFile.(Structname{1}));
%% Calling ECG Signal
% Open .mat file
[FileName,PathName] = uigetfile('*.mat','Select mat file');
if FileName==0,
return;
end
MatFile=load(fullfile(PathName,FileName));   %# pass file path as string
Structname = fieldnames(MatFile);
assignin('base', 'ecg3', MatFile.(Structname{1}));
%% Parameter entering
ecgsignal= [ecg1;ecg2;ecg3];
fs=1000;
%% Baseline Drift Removing 
% Remove the moving artifacts
ecg_signal=ecgsignal-mean(ecgsignal);
[c,l]=wavedec(ecg_signal,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
ecgcorrected=ecg_signal-a8;
%% Filtering noise
% a=fir1(20,45/100,'low');
[B,A] = butter(10,45/500,'low');
ecgfiltered_last=filter(B,A,ecgcorrected);
wt=cwt(ecgfiltered_last,50:340,'db6','scal');
f=linspace(-500,500,length(ecgfiltered_last));
figure
plot(f,fftshift(abs(fft(ecgfiltered_last))))
xlabel('freq')
ylabel('amplit')
title('frequency rep. of ecgfiltered last')
axis([-2 55 0 12e6])
F = scal2frq(50:340,'db6',1/1000);
figure
semilogy(50:340,F)
title('skala to freq')
xlabel('skala')
ylabel('frequency')
figure
plot(f,fftshift(abs(fft(ecgsignal))))
axis([-2 55 0 12e6])
% sc=wscalogram('image',wt);
% figure()
% plot(ecgfiltered_last)
toc