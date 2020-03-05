%% Calling ECG Signal
% Open .mat file
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
% %% Calling ECG Signal
% % Open .mat file
% [FileName,PathName] = uigetfile('*.mat','Select mat file');
% if FileName==0,
% return;
% end
% MatFile=load(fullfile(PathName,FileName));   %# pass file path as string
% Structname = fieldnames(MatFile);
% assignin('base', 'ecg3', MatFile.(Structname{1}));
%% Parameter entering
ecgsignal= [ecg1;ecg2];
fs=400;
bpm=100;
%% Baseline Drift Removing 
% Remove the moving artifacts
ecg_signal=ecgsignal-mean(ecgsignal);
[c,l]=wavedec(ecg_signal,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
ecgcorrected=ecg_signal-a8;
%% Filtering noise
% N=100;
% if(mod(length(ecgcorrected),2)~=0)
% w=linspace(0,fs/2,ceil(length(ecgcorrected)/2));
% wn=48;
% e=1+power(w/wn,2*N);
% for i=1:length(e)
% h(i)=abs(sqrt(1/e(i)));
% end
% h1=h(:,end:-1:1);
% h2=[h1(1:end-1) h];
% butterfilter=h2';
% ecgfft=fftshift(fft(ecgcorrected(1:end)));
% ecgfilt1=butterfilter'.*ecgfft;
% ecgfiltered=ifftshift(ecgfilt1);
% ecgfiltered_last=real(ifft(ecgfiltered));
% else
% w=linspace(0,fs/2,(length(ecgcorrected)/2));
% wn=48;
% e=1+power(w/wn,2*N);
% for i=1:length(e)
% h(i)=abs(sqrt(1/e(i)));
% end
% h1=h(:,end:-1:1);
% h2=[h1(1:end-1) h];
% butterfilt=h2';
% butterfilter=[butterfilt(1:end);zeros(1,1)];
% ecgfft=fftshift(fft(ecgcorrected(1:end)));
% ecgfilt1=butterfilter.*ecgfft;
% ecgfiltered=ifftshift(ecgfilt1);
% ecgfiltered_last=real(ifft(ecgfiltered));
% end
% ecgfiltered_last=ecgfiltered_last/(max(ecgfiltered_last));
a=fir1(100,0.192,'low');
ecgfiltered_last=filter(a,1,ecgcorrected);
wt=cwt(ecgfiltered_last,1:400,'db6','scal');
% sc=wscalogram('image',wt);
% figure()
% plot(ecgfiltered_last)
