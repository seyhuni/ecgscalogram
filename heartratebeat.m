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
bpm=400;
%% Baseline Drift Removing 
% Remove the moving artifacts
tic
ecg_signal=ecgsignal-mean(ecgsignal);
[c,l]=wavedec(ecg_signal,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
ecgcorrected=ecg_signal-a8;
%% Filtering noise
N=100;
if(mod(length(ecgcorrected),2)~=0)
w=linspace(0,fs/2,ceil(length(ecgcorrected)/2));
wn=48;
e=1+power(w/wn,2*N);
for i=1:length(e)
h(i)=abs(sqrt(1/e(i)));
end
h1=h(:,end:-1:1);
h2=[h1(1:end-1) h];
butterfilter=h2';
ecgfft=fftshift(fft(ecgcorrected(1:end)));
ecgfilt1=butterfilter'.*ecgfft;
ecgfiltered=ifftshift(ecgfilt1);
ecgfiltered_last=real(ifft(ecgfiltered));
else
w=linspace(0,fs/2,(length(ecgcorrected)/2));
wn=48;
e=1+power(w/wn,2*N);
for i=1:length(e)
h(i)=abs(sqrt(1/e(i)));
end
h1=h(:,end:-1:1);
h2=[h1(1:end-1) h];
butterfilt=h2';
butterfilter=[butterfilt(1:end);zeros(1,1)];
ecgfft=fftshift(fft(ecgcorrected(1:end)));
ecgfilt1=butterfilter.*ecgfft;
ecgfiltered=ifftshift(ecgfilt1);
ecgfiltered_last=real(ifft(ecgfiltered));
end
ecgfiltered_last=ecgfiltered_last/(max(ecgfiltered_last));

%% Wavelet Transform
% Using Wavelet Transform to Decompose signal
[c,l] = wavedec(ecgfiltered_last,8,'db6');
%% Wavelet Reconstruction with Coefficients
% DESCRIPTIVE TEXT
for t=1:8
    D(:,t)=wrcoef('d',c,l,'db6',t);
end
a2=wrcoef('a',c,l,'db6',2);
%% R-Peak Detection
% Use selected coefficients to R Peak Detection
e1=D(:,3)+D(:,4)+D(:,5);
e2=(D(:,4).*(D(:,3)+D(:,5)))/2.^8;
R_Peak_Detect_Ecg=e1.*e2;
R_Peak_Detect_Ecg_Positive = zeros(1,length(R_Peak_Detect_Ecg));
for k=1:length(R_Peak_Detect_Ecg)   
if R_Peak_Detect_Ecg(k)>0
R_Peak_Detect_Ecg_Positive(k)=R_Peak_Detect_Ecg(k);
end
end
last_ecg=R_Peak_Detect_Ecg_Positive;
threshold=max(last_ecg);
threshold=threshold*0.01;
r_peak=0;
x=ceil(fs/9);

for k=x:(length(last_ecg)-x)
        gecici=last_ecg(k-x+1:k+x);
        if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>threshold)
        r_peak(k) = last_ecg(k);
    end
end
r_peak_pos=find(r_peak>0);
for j=1:length(r_peak_pos)-1
    if(abs(r_peak_pos(j)-r_peak_pos(j+1))<=53) 
        if(r_peak(r_peak_pos(j))>r_peak(r_peak_pos(j+1)))
        r_peak(r_peak_pos(j+1))=0;
        elseif(r_peak(r_peak_pos(j))<r_peak(r_peak_pos(j+1)))
            r_peak(r_peak_pos(j))=0;
        end
    end
        
end
t=find(r_peak>0);
r_peak_pos=t;

r_peak_last=zeros(1,length(last_ecg));
for t=1:length(r_peak_pos)
    mt3=ecgfiltered_last(r_peak_pos(t)-3:r_peak_pos(t)+3);
    mn3=max(mt3);
    r_peak_last(find(mt3==mn3)+r_peak_pos(t)-4)=mn3;
end
r_peak_pos_last=find(r_peak_last>0);
for j=1:length(r_peak_pos_last)-1
    if(abs(r_peak_pos_last(j)-r_peak_pos_last(j+1))<53) 
        if(r_peak_last(r_peak_pos_last(j))>=r_peak_last(r_peak_pos_last(j+1)))
        r_peak_last(r_peak_pos_last(j+1))=0;
        elseif(r_peak_last(r_peak_pos_last(j))<=r_peak_last(r_peak_pos_last(j+1)))
            r_peak_last(r_peak_pos_last(j))=0;
        end
    end
        
end
z=find(r_peak_last>0);
r_peak_pos_last=z;
length(r_peak_pos_last)
a=zeros(1,length(ecgfiltered_last));

for i=1:length(r_peak_pos_last);
a(r_peak_pos_last(i))= ecgfiltered_last(r_peak_pos_last(i));
end
figure;
plot(ecgfiltered_last)
hold on
for i=1:length(a)
if(a(i)>0)
scatter(i,a(i),'k+')
end
end