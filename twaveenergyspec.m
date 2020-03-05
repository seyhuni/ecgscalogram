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
%% ekg sinyali iþlemleri
fs=1000;
ecgfiltered_last=[ecg1';ecg2';ecg3'];
ecgfiltered_last=ecgfiltered_last';
%% Plotting,
wt=cwt(ecgfiltered_last,1:100,'mexh','scal');
sc=wscalogram('image',wt);
en1=sum(sum(sc(:,1:60e3)))
en2=sum(sum(sc(:,60e3+1:120e3)))
en3=sum(sum(sc(:,120e3+1:180e3)))
[en1 en2 en3]
% F = scal2frq(1:100,'mexh',1/1000);
% f=linspace(-500,500,length(ecgfiltered_last));
% figure
% subplot(2,1,1)
% plot(f,fftshift(abs(fft(ecgfiltered_last))))
% axis([-2 55 0 12e3])
% xlabel('freq')
% ylabel('amplit')
% title('frequency rep. of ecgfiltered last')
% axis([-2 55 0 12e6])
% subplot(2,1,2)
% plot(f,fftshift(abs(fft(ecgsignal))))
% axis([-2 55 0 12e3])
% xlabel('freq')
% ylabel('amplt')
% title('frequency rep. of raw ecg')

% figure()
% plot(ecgfiltered_last)
% figure
% semilogy(1:100,F)
% title('scale to freq')
% xlabel('scale')
% ylabel('frequency')
% figure
toc