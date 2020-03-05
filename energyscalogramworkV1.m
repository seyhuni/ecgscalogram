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
ecgsignal=[ecg1;ecg2;ecg3];
ecgsignal=ecgsignal/(max(ecgsignal));
ecg_signal1=ecg1-mean(ecg1);
ecg_signal2=ecg2-mean(ecg2);
ecg_signal3=ecg3-mean(ecg3);
[c1,l1]=wavedec(ecg_signal1,8,'db6');
[c2,l2]=wavedec(ecg_signal2,8,'db6');
[c3,l3]=wavedec(ecg_signal3,8,'db6');
a81 = wrcoef('a',c1,l1,'db6',8);
a82 = wrcoef('a',c2,l2,'db6',8);
a83 = wrcoef('a',c3,l3,'db6',8);
ecgcorrected1=ecg_signal1-a81;
ecgcorrected2=ecg_signal2-a82;
ecgcorrected3=ecg_signal3-a83;
%% Filtering noise
N=100;
if(mod(length(ecgcorrected1),2)~=0)
w1=linspace(0,fs/2,ceil(length(ecgcorrected1)/2));
wn=46;
e1=1+power(w1/wn,2*N);
for i=1:length(e1)
h1(i)=abs(sqrt(1/e1(i)));
end
h11=h1(:,end:-1:1);
h21=[h11(1:end-1) h1];
butterfilter1=h21';
ecgfft1=fftshift(fft(ecgcorrected1(1:end)));
ecgfilt11=butterfilter1'.*ecgfft1;
ecgfiltered1=ifftshift(ecgfilt11);
ecgfiltered_last1=real(ifft(ecgfiltered1));
else
w1=linspace(0,fs/2,(length(ecgcorrected1)/2));
wn=48;
e1=1+power(w1/wn,2*N);
for i=1:length(e1)
h1(i)=abs(sqrt(1/e1(i)));
end
h11=h1(:,end:-1:1);
h21=[h1(1:end-1) h1];
butterfilt=h21';
butterfilter1=[butterfilt(1:end);zeros(1,1)];
ecgfft1=fftshift(fft(ecgcorrected1(1:end)));
ecgfilt11=butterfilter1.*ecgfft1;
ecgfiltered1=ifftshift(ecgfilt11);
ecgfiltered_last1=real(ifft(ecgfiltered1));
end
ecgfiltered_last1=ecgfiltered_last1/(max(ecgfiltered_last1));
%% Filtering noise 2
N=100;
if(mod(length(ecgcorrected2),2)~=0)
w2=linspace(0,fs/2,ceil(length(ecgcorrected2)/2));
wn=46;
e2=1+power(w2/wn,2*N);
for i=1:length(e2)
h2(i)=abs(sqrt(1/e2(i)));
end
h12=h2(:,end:-1:1);
h22=[h12(1:end-1) h2];
butterfilter2=h22';
ecgfft2=fftshift(fft(ecgcorrected2(1:end)));
ecgfilt12=butterfilter2'.*ecgfft2;
ecgfiltered2=ifftshift(ecgfilt12);
ecgfiltered_last2=real(ifft(ecgfiltered2));
else
w2=linspace(0,fs/2,(length(ecgcorrected2)/2));
wn=48;
e2=1+power(w2/wn,2*N);
for i=1:length(e2)
h2(i)=abs(sqrt(1/e2(i)));
end
h12=h2(:,end:-1:1);
h22=[h12(1:end-1) h2];
butterfilt=h22';
butterfilter2=[butterfilt(1:end);zeros(1,1)];
ecgfft2=fftshift(fft(ecgcorrected2(1:end)));
ecgfilt12=butterfilter2.*ecgfft2;
ecgfiltered2=ifftshift(ecgfilt12);
ecgfiltered_last2=real(ifft(ecgfiltered2));
end
ecgfiltered_last2=ecgfiltered_last2/(max(ecgfiltered_last2));
%% Filtering noise 3
N=100;
if(mod(length(ecgcorrected3),2)~=0)
w3=linspace(0,fs/2,ceil(length(ecgcorrected3)/2));
wn=46;
e3=1+power(w3/wn,2*N);
for i=1:length(e3)
h3(i)=abs(sqrt(1/e3(i)));
end
h13=h3(:,end:-1:1);
h23=[h13(1:end-1) h3];
butterfilter3=h23';
ecgfft3=fftshift(fft(ecgcorrected3(1:end)));
ecgfilt13=butterfilter3'.*ecgfft3;
ecgfiltered3=ifftshift(ecgfilt13);
ecgfiltered_last3=real(ifft(ecgfiltered3));
else
w3=linspace(0,fs/2,(length(ecgcorrected3)/2));
wn=48;
e3=1+power(w3/wn,2*N);
for i=1:length(e3)
h3(i)=abs(sqrt(1/e3(i)));
end	
h13=h3(:,end:-1:1);
h23=[h13(1:end-1) h3];
butterfilt=h23';
butterfilter3=[butterfilt(1:end);zeros(1,1)];
ecgfft3=fftshift(fft(ecgcorrected3(1:end)));
ecgfilt13=butterfilter3.*ecgfft3;
ecgfiltered3=ifftshift(ecgfilt13);
ecgfiltered_last3=real(ifft(ecgfiltered3));
end
ecgfiltered_last3=ecgfiltered_last3/(max(ecgfiltered_last3));
%% Final Signal
ecgfiltered_last= [ecgfiltered_last1;ecgfiltered_last2;ecgfiltered_last3];
%% Plotting,
wt=cwt(ecgfiltered_last,50:340,'db6','scal');
% F = scal2frq(50:340,'db6',1/1000);
f=linspace(-500,500,length(ecgfiltered_last));
figure
subplot(2,1,1)
plot(f,fftshift(abs(fft(ecgfiltered_last))))
axis([-2 55 0 12e3])
xlabel('freq')
ylabel('amplit')
title('frequency rep. of ecgfiltered last')
% axis([-2 55 0 12e6])
subplot(2,1,2)
plot(f,fftshift(abs(fft(ecgsignal))))
axis([-2 55 0 12e3])
xlabel('freq')
ylabel('amplt')
title('frequency rep. of raw ecg')
figure
sc=wscalogram('image',wt);
% figure()
% plot(ecgfiltered_last)
% figure
% semilogy(50:340,F)
% title('scale to freq')
% xlabel('scale')
% ylabel('frequency')
% figure
toc
en1=sum(sum(sc(:,1:60e3)))
en2=sum(sum(sc(:,60e3+1:120e3)))
en3=sum(sum(sc(:,120e3+1:180e3)))