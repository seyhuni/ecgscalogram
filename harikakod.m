%% Calling ECG Signal
% Open .mat file
srcFiles = dir('C:\Users\huseyinyanik\Desktop\akademik\EKG ÇALIÞMALARI\antrene work\Erkek Kontrol\*.mat');
for kkk=52:length(srcFiles)
% [FileName,PathName] = uigetfile('*.mat','Select mat file');
FileName=srcFiles(kkk).name;
PathName='C:\Users\huseyinyanik\Desktop\akademik\EKG ÇALIÞMALARI\antrene work\Erkek Kontrol';
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
%% Adaptive Filtering
fl=30;  % lower and upper cutoff freqs
wcl=2*fl/fs;
% faxes=linspace(-fs/2,fs/2,length(ecgsignal));
% plot(faxes,fftshift(abs(fft(ecgsignal))));
% axis([0 100 0 4e6])
N=20;
wn=[wcl];
b=fir1(N,wn,hamming(N+1)); %notch filter with hamming window
% [h,w]= freqz(b,1,256);
% j=1:length(ecgcorrected);
iv=zeros(1,N); %initialization vector for all taps to zero
ecgfiltered_last=filter(b,1,ecgcorrected,iv);
ecgfiltered_last = sgolayfilt(ecgfiltered_last,3,11);
% figure; plot(ecgsignal/max(ecgsignal)); title('Raw signal');
% figure; plot(ecgfiltered_last/max(ecgfiltered_last)); title('Filtered signal');
ecgfiltered_last=ecgfiltered_last/max(ecgfiltered_last);
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
    mt3=ecgfiltered_last(r_peak_pos(t)-60:r_peak_pos(t)+60);
    mn3=max(mt3);
    r_peak_last(find(mt3==mn3)+r_peak_pos(t)-61)=mn3;
end
r_peak_pos_last=find(r_peak_last>0);

for i=1:length(r_peak_pos_last)-1
if((r_peak_pos_last(i+1)-r_peak_pos_last(i))<500) & (ecgfiltered_last(r_peak_pos_last(i+1))>ecgfiltered_last(r_peak_pos_last(i)))
        r_peak_pos_last(i)= [];
    else if(((r_peak_pos_last(i+1)-r_peak_pos_last(i))<500) & (ecgfiltered_last(r_peak_pos_last(i+1))<ecgfiltered_last(r_peak_pos_last(i))))
                r_peak_pos_last(i+1)= [];
        end
    end
end
%% Q & S Detection
hv=ones(1,length(ecgfiltered_last));
s_peak=sqrt(-1)*hv;
q_peak=s_peak;

for k=1:length(r_peak_pos_last)
    mt=ecgfiltered_last(r_peak_pos_last(k):r_peak_pos_last(k)+40);
    mn=min(mt);
    s_peak(find(mt==mn)+r_peak_pos_last(k)-1)=mn;
end
s_peak_pos_last=find(s_peak~=sqrt(-1));

for m=1:length(r_peak_pos_last)
    mt2=ecgfiltered_last(r_peak_pos_last(m)-50:r_peak_pos_last(m));
    mn2=min(mt2);
    q_peak(r_peak_pos_last(m)-50-1+find(mt2==mn2))= mn2;
end
q_peak_pos_last=find(q_peak~=sqrt(-1));

%% P and T Detection   ******* Burda kaldým Buraya kadar Süper 
e4=D(:,4)+D(:,5)+D(:,6)+D(:,7)+D(:,8);

t_peak=hv*sqrt(-1);
for t=1:length(s_peak_pos_last)-1
    mt6=e4(s_peak_pos_last(t):s_peak_pos_last(t)+250);
    mn6=max(mt6);
    t_peak(find(mt6==mn6)+s_peak_pos_last(t)-1)=mn6;
end
t_peak_pos=find(t_peak~=sqrt(-1));

t_peak_last=hv*sqrt(-1);

for t=1:length(t_peak_pos)
    mt9=ecgfiltered_last(t_peak_pos(t)-10:t_peak_pos(t)+10);
    mn9=max(mt9);
    in = find(mt9==mn9);
    t_peak_last(in(1)+t_peak_pos(t)-11)=mn9;
end

t_peak_pos_last = find(t_peak_last~=sqrt(-1));

% for x=1:length(t_peak_pos_last)
% if(r_peak_pos_last(x+1)-t_peak_pos_last(x)<(ceil(fs/bpm)*6))
%     t_peak_last(t_peak_pos_last(x))=0;
% mt13=ecgfiltered_last(s_peak_pos_last(x):s_peak_pos_last(x)+(ceil(fs/bpm)*8));
% mn13=max(mt13);
% in=find(mt13==mn13);
% t_peak_last(in(1)+s_peak_pos_last(x)-1)=mn13;
% end
% end
% t_peak_pos_last = find(t_peak_last~=sqrt(-1));

p_peak=hv*sqrt(-1);

for p=2:length(q_peak_pos_last)
    mt10=e4(q_peak_pos_last(p)- 120:q_peak_pos_last(p));
    mn10=max(mt10);
    p_peak(find(mt10==mn10)+q_peak_pos_last(p)-120-1)=mn10;
end

p_peak_pos=find(p_peak~=sqrt(-1));

p_peak_last=hv*sqrt(-1);

for p=1:length(p_peak_pos)
    mt11=ecgfiltered_last(p_peak_pos(p)- 10:p_peak_pos(p)+10);
    mn11=max(mt11);
    in = find(mt11==mn11);
    p_peak_last(in(1)+p_peak_pos(p)-10+1)=mn11;
end

p_peak_pos_last=find(p_peak_last~=sqrt(-1));
%%
%%P ve T Start and Final Points
% P Start
p_peak_start=hv*sqrt(-1);
for i=2:length(p_peak_pos_last)
    mt15=ecgfiltered_last(p_peak_pos_last(i)-80:p_peak_pos_last(i));
    mn15=min(mt15);
    in = find(mt15==mn15);
    p_peak_start(in(1)+p_peak_pos_last(i)-80-1)=mn15;
end
p_peak_start_pos=find(p_peak_start~=sqrt(-1));

%P Final
p_peak_final=hv*sqrt(-1);
for i=1:length(p_peak_pos_last)
    mt16=ecgfiltered_last(p_peak_pos_last(i):p_peak_pos_last(i)+45);
    mn16=min(mt16);
    in = find(mt16==mn16);
    p_peak_final(in(1)+p_peak_pos_last(i)-1)=mn16;
end
p_peak_final_pos=find(p_peak_final~=sqrt(-1));

for i=1:length(p_peak_pos_last)-1
    if(q_peak_pos_last(i+1)<=p_peak_final_pos(i))
        q_peak_pos_last(i+1)=q_peak_pos_last(i+1)+3;
        p_peak_final_pos(i)=p_peak_final_pos(i)-3;
    end
end

% T Final
t_peak_final=hv*sqrt(-1);
for i=1:length(t_peak_pos_last)
    mt18=ecgfiltered_last(t_peak_pos_last(i):t_peak_pos_last(i)+110);
    mn18=min(mt18);
    in = find(mt18==mn18);
    t_peak_final(in(1)+t_peak_pos_last(i)-1)=mn18;
end
t_peak_final_pos=find(t_peak_final~=sqrt(-1));

    t_peak_start=hv*sqrt(-1);   
for i=2:length(t_peak_pos_last)
    mt17=ecgfiltered_last(t_peak_pos_last(i)-120:t_peak_pos_last(i));
    mn17=min(mt17);
    in = find(mt17==mn17);
    t_peak_start(in(1)+t_peak_pos_last(i)-121)=mn17;
end
t_peak_start_pos=find(t_peak_start~=sqrt(-1));     
%% BPM Calculating
% Calculate Beats per minute
% Peakleri Ustune Bindirme
% a=zeros(1,length(ecgfiltered_last));
% b=hv*sqrt(-1);
% c=b;
% d=b;
% e=b;
% ps=b;
% pf=b;
% ts=b;
% tf=b;
% for i=1:length(r_peak_pos_last);
% a(r_peak_pos_last(i))= ecgfiltered_last(r_peak_pos_last(i));
% end
% for i=1:length(s_peak_pos_last);
% b(s_peak_pos_last(i))= ecgfiltered_last(s_peak_pos_last(i));
% end
% for i=1:length(q_peak_pos_last);
% c(q_peak_pos_last(i))= ecgfiltered_last(q_peak_pos_last(i));
% end
% for i=1:length(t_peak_pos_last);
% d(t_peak_pos_last(i))=ecgfiltered_last(t_peak_pos_last(i));
% end
% for i=1:length(p_peak_pos_last);
% e(p_peak_pos_last(i))= ecgfiltered_last(p_peak_pos_last(i));
% end
% for i=1:length(p_peak_start_pos);
% ps(p_peak_start_pos(i))=ecgfiltered_last(p_peak_start_pos(i));
% end
% for i=1:length(p_peak_final_pos);
% pf(p_peak_final_pos(i))=ecgfiltered_last(p_peak_final_pos(i));
% end
% for i=1:length(t_peak_start_pos);
% ts(t_peak_start_pos(i))=ecgfiltered_last(t_peak_start_pos(i));
% end
% for i=1:length(t_peak_final_pos);
% tf(t_peak_final_pos(i))=ecgfiltered_last(t_peak_final_pos(i));
% end
figure;
hold on;
plot(ecgfiltered_last)
plot(r_peak_pos_last,ecgfiltered_last(r_peak_pos_last),'r+','MarkerFaceColor','r')
plot(s_peak_pos_last,ecgfiltered_last(s_peak_pos_last),'r*','MarkerFaceColor','r')
plot(q_peak_pos_last,ecgfiltered_last(q_peak_pos_last),'r.','MarkerFaceColor','r')
plot(p_peak_start_pos,ecgfiltered_last(p_peak_start_pos),'k+','MarkerFaceColor','r')
plot(p_peak_pos_last,ecgfiltered_last(p_peak_pos_last),'k*','MarkerFaceColor','r')
plot(p_peak_final_pos,ecgfiltered_last(p_peak_final_pos),'k.','MarkerFaceColor','r')
plot(t_peak_start_pos,ecgfiltered_last(t_peak_start_pos),'b+','MarkerFaceColor','r')
plot(t_peak_pos_last,ecgfiltered_last(t_peak_pos_last),'b*','MarkerFaceColor','r')
plot(t_peak_final_pos,ecgfiltered_last(t_peak_final_pos),'b.','MarkerFaceColor','r')
pause(10)

% plot(ecgfiltered_last)
% hold on
% for i=1:length(a)
% if(a(i)>0)
% scatter(i,a(i),'k+')
% end
% end
% for i=1:length(b)
% if(b(i)~=sqrt(-1))
% scatter(i,b(i),'r+')
% end
% end
% for i=1:length(c)
% if(c(i)~=sqrt(-1))
% scatter(i,c(i),'b+')
% end
% end
% for i=1:length(d)
% if(d(i)~=sqrt(-1))
% scatter(i,d(i),'g+')
% end
% end
% for i=1:length(e)
% if(e(i)~=sqrt(-1))
% scatter(i,e(i),'r*')
% end
% end
% for i=1:length(ps)
% if(ps(i)~=sqrt(-1))
% scatter(i,ps(i),'k*')
% end
% end
% for i=1:length(pf)
% if(pf(i)~=sqrt(-1))
% scatter(i,pf(i),'g*')
% end
% end
% for i=1:length(ts)
% if(ts(i)~=sqrt(-1))
% scatter(i,ts(i),'ms')
% end
% end
% for i=1:length(tf)
% if(tf(i)~=sqrt(-1))
% scatter(i,tf(i),'m+')
% end
% end
%% ST Interval and T Peak 

st_interval= zeros(1,length(t_peak_start_pos));
for i =1:length(st_interval)
    st_interval(i)=t_peak_final_pos(i)-s_peak_pos_last(i); 
    %t bitiþ - s
end
mean_st_interval= mean(st_interval)

st_segment=zeros(1,length(t_peak_pos_last));
for i =1:length(st_interval)
    st_segment(i)=t_peak_start_pos(i)-s_peak_pos_last(i);
    % t baþlangýç - s
end
mean_st_segment=mean(st_segment);
t_peak_amplitude=zeros(1,length(t_peak_pos_last));
for i=1:length(t_peak_amplitude)-1
    t_peak_amplitude(i) = ecgfiltered_last(t_peak_pos_last(i))-ecgfiltered_last(t_peak_final_pos(i));
end
t_peak_amp_mean=mean(t_peak_amplitude)

%% QRS Interval and amplitude
qrs_amp=zeros(1,length(r_peak_pos_last));
qrs_int=zeros(1,length(q_peak_pos_last));
for i=1:length(qrs_int)
    qrs_int(i)=s_peak_pos_last(i)-q_peak_pos_last(i);
end
for i=1:length(qrs_amp)
    qrs_amp(i)=ecgfiltered_last(r_peak_pos_last(i))-ecgfiltered_last(s_peak_pos_last(i));
end
mean_qrs_int=mean(qrs_int)
mean_qrs_amp=mean(qrs_amp)
%% RR interval
rr_int=zeros(1,length(r_peak_pos_last)-1);

for i=1:length(r_peak_pos_last)-1
    rr_int(i)=(r_peak_pos_last(i+1)-r_peak_pos_last(i));
end
mean_rr=mean(rr_int)
%% QT interval 
qt_int=zeros(1,length(q_peak_pos_last));
if(t_peak_final_pos(1)>q_peak_pos_last(1))
for i=1:length(q_peak_pos_last)-2
    qt_int(i)=t_peak_final_pos(i)-q_peak_pos_last(i);
end
else
    for i=1:length(q_peak_post_last)-1
        qt_int(i)=t_peak_final_pos(i+1)-q_peak_pos_last(i);
    end
end
qt_intmean=mean(qt_int)

%% PR interval and P amp
p_amp=zeros(1,length(p_peak_pos_last));
pr_int=zeros(1,length(p_peak_pos_last));
for i=1:length(p_peak_pos_last)-1
    p_amp(i)=ecgfiltered_last(p_peak_pos_last(i))-ecgfiltered_last(p_peak_start_pos(i));
end
p_ampmean=mean(p_amp)

if(q_peak_pos_last(1)>p_peak_start_pos(1))
for i=1:length(p_peak_pos_last)-1
    pr_int(i)=q_peak_pos_last(i)-p_peak_start_pos(i);
end
else if(q_peak_pos_last(2)>p_peak_start_pos(1))
    for i=1:length(p_peak_pos_last)-1
        pr_int(i)=q_peak_pos_last(i+1)-p_peak_start_pos(i);
    end
    else
        for i=1:length(p_peak_pos_last)-1
        pr_int(i)=q_peak_pos_last(i+2)-p_peak_start_pos(i);
    end
    end
end
pr_mean=mean(pr_int)

basefilename=srcFiles(kkk).name;
dlmwrite([basefilename '.csv'] ,[length(r_peak_pos_last)' mean_rr' mean_qrs_int' mean_qrs_amp' t_peak_amp_mean' mean_st_segment' mean_st_interval' qt_intmean' p_ampmean' pr_mean'],'precision','%20.5f');
end
%% RR Ýnterval
% KARAKTERÝSTK SÜRE VE GENLÝK DEÐERLERÝ
% %% RR Interval
% 
% toplam=0;
% 
% for i=1:length(r_peak_pos_last)-1
%     toplam=toplam+(r_peak_pos_last(i+1)-r_peak_pos_last(i));
% end
% 
% rraverage=toplam/(length(r_peak_pos_last)-1);
% bradikardlim=rraverage/80*60;
% tasikardlim=rraverage/80*100;
% place=zeros(1,length(r_peak_pos_last)-1);
% for i=1:length(r_peak_pos_last)-1
%     if(((r_peak_pos_last(i+1)-r_peak_pos_last(i))<bradikardlim) | ((r_peak_pos_last(i+1)-r_peak_pos_last(i))>tasikardlim))
%       
%         place(i)=i;
%         
% %     elseif(((r_peak_pos_last(i+1)-r_peak_pos_last(i))<(bradikardlim+32) | ((r_peak_pos_last(i+1)-r_peak_pos_last(i))>(tasikardlim-32))))
% %        
% %         place(i)=i;
%         
%     else
%         
%         place(i)=0;
%         
%     end
% end
% 
% arrythmi = find(place~=0);
% if(arrythmi~=0)
% for i=1:length(arrythmi)
%     display(sprintf('heart has abnormalities because of %d. R peak',arrythmi(i)))
% end
% else
%     display('heart has no arrythmia');
% end
% 
% %% Veri kayit islemi
% % R Kayit 
% 
% r_peak_pos_sec=zeros(1,length(r_peak_pos_last));
% 
% for i=1:length(r_peak_pos_last)
%     r_peak_pos_sec(i)=r_peak_pos_last(i)/fs;
% end
% 
% genlikr = zeros(1,length(r_peak_pos_last));
% 
% for i=1:length(r_peak_pos_last)
%     genlikr(i)=ecgsignal(r_peak_pos_last(i));
% end
% 
% dlmwrite('Rbilgi.csv',[r_peak_pos_last' r_peak_pos_sec' genlikr'],'precision','%20.5f');
% 
% % Q Kayýt
% 
% q_peak_pos_sec=zeros(1,length(q_peak_pos_last));
% 
% for i=1:length(q_peak_pos_last)
%     q_peak_pos_sec(i)=q_peak_pos_last(i)/fs;
% end
% 
% genlikq = zeros(1,length(q_peak_pos_last));
% 
% for i=1:length(q_peak_pos_last)
%     genlikq(i)=ecgsignal(q_peak_pos_last(i));
% end
% 
% dlmwrite('Qbilgi.csv',[q_peak_pos_last' q_peak_pos_sec' genlikq'],'precision','%20.5f');
% 
% % S Kayýt
% 
% s_peak_pos_sec=zeros(1,length(s_peak_pos_last));
% 
% for i=1:length(s_peak_pos_last)
%     s_peak_pos_sec(i)=s_peak_pos_last(i)/fs;
% end
% 
% genliks = zeros(1,length(s_peak_pos_last));
% 
% for i=1:length(s_peak_pos_last)
%     genliks(i)=ecgsignal(s_peak_pos_last(i));
% end
% 
% dlmwrite('Sbilgi.csv',[s_peak_pos_last' s_peak_pos_sec' genliks'],'precision','%20.5f');
% 
% %P Kayýt
% 
% p_peak_pos_sec=zeros(1,length(p_peak_pos_last));
% 
% for i=1:length(p_peak_pos_last)
%     p_peak_pos_sec(i)=p_peak_pos_last(i)/fs;
% end
% 
% genlikp = zeros(1,length(p_peak_pos_last));
% 
% for i=1:length(p_peak_pos_last)
%     genlikp(i)=ecgsignal(p_peak_pos_last(i));
% end
% 
% dlmwrite('Pbilgi.csv',[p_peak_pos_last' p_peak_pos_sec' genlikp'],'precision','%20.5f');
% 
% % T Kayýt
% 
% t_peak_pos_sec=zeros(1,length(t_peak_pos_last));
% 
% for i=1:length(t_peak_pos_last)
%     t_peak_pos_sec(i)=t_peak_pos_last(i)/fs;
% end
% 
% genlikt = zeros(1,length(t_peak_pos_last));
% 
% for i=1:length(t_peak_pos_last)
%     genlikt(i)=ecgsignal(t_peak_pos_last(i));
% end
% 
% dlmwrite('Tbilgi.csv',[t_peak_pos_last' t_peak_pos_sec' genlikt'],'precision','%20.5f');