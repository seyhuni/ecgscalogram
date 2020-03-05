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
    if(abs(r_peak_pos(j)-r_peak_pos(j+1))<=23) 
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
    if(abs(r_peak_pos_last(j)-r_peak_pos_last(j+1))<23) 
        if(r_peak_last(r_peak_pos_last(j))>=r_peak_last(r_peak_pos_last(j+1)))
        r_peak_last(r_peak_pos_last(j+1))=0;
        elseif(r_peak_last(r_peak_pos_last(j))<=r_peak_last(r_peak_pos_last(j+1)))
            r_peak_last(r_peak_pos_last(j))=0;
        end
    end
        
end
z=find(r_peak_last>0);
r_peak_pos_last=z;
%% Q & S Detection
e3=D(:,2)+D(:,3)+D(:,4)+D(:,5);
for i=1:length(e3)
    e(i+3)=e3(i);
end
e=[zeros(1,4) e zeros(1,4)];
for i=3:length(e3)
    f(i)=(-e(i+2)+8*(e(i+1))-8*(e(i-1))+e(i-2))/12;
end
S_Peak_Detect_Ecg_Negative = zeros(1,length(f));
for k=1:length(f)   
if f(k)<0
S_Peak_Detect_Ecg_Negative(k)=f(k);
end
end
last_ecg_s=S_Peak_Detect_Ecg_Negative / max(abs(S_Peak_Detect_Ecg_Negative));
last_ecg_s=abs(last_ecg_s);
threshold_s=ones(1,length(last_ecg_s));
threshold_s=threshold_s*0.175;
hv=ones(1,length(last_ecg_s));
s_peak=sqrt(-1)*hv;
q_peak=s_peak;
sbeatcount=0;

for k=1:length(r_peak_pos_last)
    mt=a2(r_peak_pos_last(k):r_peak_pos_last(k)+ceil(fs/bpm)*2);
    mn=min(mt);
    s_peak(find(mt==mn)+r_peak_pos_last(k)-1)=mn;
    sbeatcount=sbeatcount+1;
end
s_peak_pos=find(s_peak~=sqrt(-1));

for m=1:length(r_peak_pos_last)
    mt2=a2(r_peak_pos_last(m)-(ceil(fs/bpm)*2):r_peak_pos_last(m));
    mn2=min(mt2);
    q_peak(r_peak_pos_last(m)-(ceil(fs/bpm)*2)-1+find(mt2==mn2))= mn2;
end
q_peak_pos=find(q_peak~=sqrt(-1));
s_peak_last=hv*sqrt(-1);
for t=1:length(s_peak_pos)
    mt7=ecgfiltered_last(s_peak_pos(t)-(ceil(fs/bpm)+3):s_peak_pos(t)+(ceil(fs/bpm)+3));
    mn7=min(mt7);
    in = find(mt7==mn7);
    s_peak_last(in(1)+s_peak_pos(t)-(ceil(fs/bpm)+4))=mn7;
%     genlik = s_peak_last(find(mt7==mn7)+s_peak_pos(t)-(ceil(fs/bpm)+4))
%     index = find(mt7==mn7)+s_peak_pos(t)-(ceil(fs/bpm)+4)
    
%     pause
    
end
s_peak_pos_last=find(s_peak_last~=sqrt(-1));
q_peak_last=hv*sqrt(-1);
for t=1:length(q_peak_pos)
    mt8=ecgfiltered_last(q_peak_pos(t)-(ceil(fs/bpm)+1):q_peak_pos(t)+(ceil(fs/bpm)+1));
    mn8=min(mt8);
    in = find(mt8==mn8);
    q_peak_last(in(1)+q_peak_pos(t)-(ceil(fs/bpm)+2))=mn8;
end
q_peak_pos_last=find(q_peak_last~=sqrt(-1));
%% P and T Detection   ******* Burda kaldým Buraya kadar Süper 
e4=D(:,4)+D(:,5)+D(:,6)+D(:,7)+D(:,8);

t_peak=hv*sqrt(-1);
for t=1:length(s_peak_pos_last)-1
    mt6=e4(s_peak_pos_last(t)+ceil(fs/bpm)*2:s_peak_pos_last(t)+ceil(fs/bpm)*30);
    mn6=max(mt6);
    t_peak(find(mt6==mn6)+s_peak_pos_last(t)+ceil(fs/bpm)*2-1)=mn6;
end
t_peak_pos=find(t_peak~=sqrt(-1));

t_peak_last=hv*sqrt(-1);

for t=1:length(t_peak_pos)
    mt9=ecgfiltered_last(t_peak_pos(t)-2:t_peak_pos(t)+ceil(fs/bpm));
    mn9=max(mt9);
    in = find(mt9==mn9);
    t_peak_last(in(1)+t_peak_pos(t)-3)=mn9;
end

t_peak_pos_last = find(t_peak_last~=sqrt(-1));

for x=1:length(t_peak_pos_last)
if(r_peak_pos_last(x+1)-t_peak_pos_last(x)<(ceil(fs/bpm)*3))
    t_peak_last(t_peak_pos_last(x))=0;
mt13=ecgfiltered_last(s_peak_pos_last(x):s_peak_pos_last(x)+(ceil(fs/bpm)*8));
mn13=max(mt13);
in=find(mt13==mn13);
t_peak_last(in(1)+s_peak_pos_last(x)-1)=mn13;
end
end
t_peak_pos_last = find(t_peak_last~=sqrt(-1));

p_peak=hv*sqrt(-1);

for p=2:length(q_peak_pos_last)
    mt10=e4(q_peak_pos_last(p)- ceil((fs/bpm)*9):q_peak_pos_last(p));
    mn10=max(mt10);
    p_peak(find(mt10==mn10)+q_peak_pos_last(p)-ceil((fs/bpm)*9)-1)=mn10;
end

p_peak_pos=find(p_peak~=sqrt(-1));

p_peak_last=hv*sqrt(-1);

for p=2:length(p_peak_pos)
    mt11=ecgfiltered_last(p_peak_pos(p)- ceil((fs/bpm)*10):p_peak_pos(p)+ceil(fs/bpm)-1);
    mn11=max(mt11);
    in = find(mt11==mn11);
    p_peak_last(in(1)+p_peak_pos(p)-(ceil((fs/bpm)*10)+1))=mn11;
end

p_peak_pos_last=find(p_peak_last~=sqrt(-1));

for x=1:length(p_peak_pos_last)
if(p_peak_pos_last(x)-q_peak_pos_last(x+1)>=0)
    p_peak_last(p_peak_pos_last(x))=sqrt(-1);
   mt12=ecgfiltered_last(q_peak_pos_last(x+1)-(ceil(fs/bpm)*8):q_peak_pos_last(x+1));
   mn12=max(mt12);
   in= find(mt12==mn12);
   p_peak_last(in(1)+q_peak_pos_last(x+1)-(ceil(fs/bpm)*8+1))=mn12;
end
end
p_peak_pos_last=find(p_peak_last~=sqrt(-1));

for x=1:length(p_peak_pos_last)
if(q_peak_pos_last(x+1)-p_peak_pos_last(x)<=ceil(fs/bpm)+9)
    p_peak_last(p_peak_pos_last(x))=sqrt(-1);
   mt14=ecgfiltered_last(q_peak_pos_last(x+1)-(ceil(fs/bpm)*7):q_peak_pos_last(x+1)-(ceil(fs/bpm)*5));
   mn14=max(mt14);
   in= find(mt14==mn14);
   p_peak_last(in(1)+q_peak_pos_last(x+1)-(ceil(fs/bpm)*7+1))=mn14;
end
end
p_peak_pos_last=find(p_peak_last~=sqrt(-1));

%%
%%P ve T Start and Final Points
% P Start
p_peak_start=hv*sqrt(-1);
for i=1:length(p_peak_pos_last)
    mt15=ecgfiltered_last(p_peak_pos_last(i)-ceil(fs/bpm)*4:p_peak_pos_last(i));
    mn15=min(mt15);
    in = find(mt15==mn15);
    p_peak_start(in(1)+p_peak_pos_last(i)-ceil(fs/bpm)*4-1)=mn15;
end
p_peak_start_pos=find(p_peak_start~=sqrt(-1));

%P Final
p_peak_final=hv*sqrt(-1);
for i=1:length(p_peak_pos_last)
    mt16=ecgfiltered_last(p_peak_pos_last(i):p_peak_pos_last(i)+ceil(fs/bpm)*2);
    mn16=min(mt16);
    in = find(mt16==mn16);
    p_peak_final(in(1)+p_peak_pos_last(i)-1)=mn16;
end
p_peak_final_pos=find(p_peak_final~=sqrt(-1));

for i=1:length(p_peak_pos_last)
    if(q_peak_pos_last(i+1)<=p_peak_final_pos(i))
        q_peak_pos_last(i+1)=q_peak_pos_last(i+1)+4;
        p_peak_final_pos(i)=p_peak_final_pos(i)-4;
    end
end

% % % % % % % T Start
% % % % % % t_peak_start=hv*sqrt(-1);
% % % % % % for i=1:length(t_peak_pos_last)
% % % % % %     mt17=ecgsignal(t_peak_pos_last(i)-ceil(fs/bpm)*9:t_peak_pos_last(i));
% % % % % %     mn17=min(mt17);
% % % % % %     in = find(mt17==mn17);
% % % % % %     t_peak_start(in(1)+t_peak_pos_last(i)-ceil(fs/bpm)*9-1)=mn17;
% % % % % % end
% % % % % % t_peak_start_pos=find(t_peak_start~=sqrt(-1));

% for i=1:length(t_peak_start_pos)
%     if(t_peak_start_pos(i)<=s_peak_pos_last(i))
%         t_peak_start_pos(i)=s_peak_pos_last(i)+round(fs/bpm)*2-1;
%     end
% end

% T Final
t_peak_final=hv*sqrt(-1);
for i=1:length(t_peak_pos_last)
    mt18=ecgfiltered_last(t_peak_pos_last(i):t_peak_pos_last(i)+28);
    mn18=min(mt18);
    in = find(mt18==mn18);
    t_peak_final(in(1)+t_peak_pos_last(i)-1)=mn18;
end
t_peak_final_pos=find(t_peak_final~=sqrt(-1));

% for i=1:length(t_peak_final_pos)-2
%     if(t_peak_final_pos(i)>=p_peak_start_pos(i))
%         t_peak_final_pos(i)=t_peak_final_pos(i)-(t_peak_final_pos(i)-p_peak_start_pos(i)+2);
%         p_peak_start_pos(i)=p_peak_start_pos(i)+(t_peak_final_pos(i)-p_peak_start_pos(i)+2);
%     end
% end
    t_peak_start=hv*sqrt(-1);   %%%%%%%%%%%%%%% Burda kaldým t_peak start ayarlamakta sýkýntý yaþýyorum
    %%%%%%%%%%%%%%%%% Çünkü ayný nokta yok. Yakýn noktayý bulmak
    %%%%%%%%%%%%%%%%% zorundayým.
for i=1:length(t_peak_pos_last)
    mt17=ecgfiltered_last(t_peak_pos_last(i)-ceil(fs/bpm)*4:t_peak_pos_last(i));
    mn17=min(mt17);
    in = find(mt17==mn17);
    t_peak_start(in(1)+t_peak_pos_last(i)-ceil(fs/bpm)*4-1)=mn17;
end
t_peak_start_pos=find(t_peak_start~=sqrt(-1));
% for i=1:length(s_peak_pos_last)-1
%     if(ecgfiltered_last(s_peak_pos_last(i))<ecgfiltered_last(t_peak_final_pos(i)))
%         for i=1:length(t_peak_pos_last)
%             mt17=ecgfiltered_last(t_peak_pos_last(i)-12:t_peak_pos_last(i));
% %benzerlik=ecgfiltered_last(t_peak_final_pos(i))-40:0.0001:ecgfiltered_last(t_peak_final_pos(i)):0.0001:ecgfiltered_last(t_peak_final_pos(i)+40);
%             val=ecgfiltered_last(t_peak_final_pos(i));
%             temp=abs(mt17-val);
%             [idx idx] = min(temp);
%             closest=mt17(idx);
%             t_peak_start(i)=closest;
%             t_peak_start_pos(i)=find(t_peak_start(i)==ecgfiltered_last);
%         end
%         
% %         benzerlik=ecgfiltered_last(t_peak_final_pos(i))-40:1:ecgfiltered_last(t_peak_final_pos(i)):1:ecgfiltered_last(t_peak_final_pos(i)+8);
% %         benzerlik=benzerlik';
%         %%%%find(mt17(17)==benzerlik(1:length(benzerlik)))                
%     end
%     
% %     benzerlik=ecgfiltered_last(t_peak_final_pos(2))-40:0.01:ecgfiltered_last(t_peak_final_pos(2)):0.01:ecgfiltered_last(t_peak_final_pos(2)+40);
%     if(ecgfiltered_last(s_peak_pos_last(i))>ecgfiltered_last(t_peak_final_pos(i)))
%     for i=1:length(t_peak_pos_last)
%         mt17=ecgfiltered_last(t_peak_pos_last(i):t_peak_pos_last(i)+18);
%         val=ecgfiltered_last(s_peak_pos_last(i));
%         temp=abs(mt17-val);
%         [idx idx] = min(temp);
%         closest=mt17(idx);
%         s_peak_last(i)=closest;
%         t_peak_final_pos(i)=find(s_peak_last(i)==ecgfiltered_last);
%     end
%      end
% end
        
        
%% BPM Calculating
% Calculate Beats per minute
% Peakleri Ustune Bindirme
a=zeros(1,length(ecgfiltered_last));
b=hv*sqrt(-1);
c=b;
d=b;
e=b;
ps=b;
pf=b;
ts=b;
tf=b;
for i=1:length(r_peak_pos_last);
a(r_peak_pos_last(i))= ecgfiltered_last(r_peak_pos_last(i));
end
for i=1:length(s_peak_pos_last);
b(s_peak_pos_last(i))= ecgfiltered_last(s_peak_pos_last(i));
end
for i=1:length(q_peak_pos_last);
c(q_peak_pos_last(i))= ecgfiltered_last(q_peak_pos_last(i));
end
for i=1:length(t_peak_pos_last);
d(t_peak_pos_last(i))=ecgfiltered_last(t_peak_pos_last(i));
end
for i=1:length(p_peak_pos_last);
e(p_peak_pos_last(i))= ecgfiltered_last(p_peak_pos_last(i));
end
for i=1:length(p_peak_start_pos);
ps(p_peak_start_pos(i))=ecgfiltered_last(p_peak_start_pos(i));
end
for i=1:length(p_peak_final_pos);
pf(p_peak_final_pos(i))=ecgfiltered_last(p_peak_final_pos(i));
end
for i=1:length(t_peak_start_pos);
ts(t_peak_start_pos(i))=ecgfiltered_last(t_peak_start_pos(i));
end
for i=1:length(t_peak_final_pos);
tf(t_peak_final_pos(i))=ecgfiltered_last(t_peak_final_pos(i));
end
figure;
plot(ecgfiltered_last)
hold on
for i=1:length(a)
if(a(i)>0)
scatter(i,a(i),'k+')
end
end
for i=1:length(b)
if(b(i)~=sqrt(-1))
scatter(i,b(i),'ko')
end
end
for i=1:length(c)
if(c(i)~=sqrt(-1))
scatter(i,c(i),'bo')
end
end
for i=1:length(d)
if(d(i)~=sqrt(-1))
scatter(i,d(i),'r+')
end
end
for i=1:length(e)
if(e(i)~=sqrt(-1))
scatter(i,e(i),'r*')
end
end
for i=1:length(ps)
if(ps(i)~=sqrt(-1))
scatter(i,ps(i),'k*')
end
end
for i=1:length(pf)
if(pf(i)~=sqrt(-1))
scatter(i,pf(i),'g*')
end
end
for i=1:length(ts)
if(ts(i)~=sqrt(-1))
scatter(i,ts(i),'ms')
end
end
for i=1:length(tf)
if(tf(i)~=sqrt(-1))
scatter(i,tf(i),'gd')
end
end
% T Peak FFT leri
% if(length(t_peak_start_pos)>250) 
% ecgt=0;
% ecgt1=0;
% % for i=1:8
% %     ecgt=ecgsignal(t_peak_start_pos(i):t_peak_final_pos(i));
% %     ecgt1(1:length(ecgt),i)=ecgt;
% %     figure;
% %     plot(fftshift(abs(fft(ecgt1(:,i)))));
% % end
% 
% % Burda ECG sinyalinin fft lerin içizdiriyorum. 10-40-70 etc.
% % sinyal çizdirme için
% 
% for i=0:9
%     ecgt=ecgfiltered_last(t_peak_start_pos(i*30+10):t_peak_final_pos(i*30+10));
%     ecgt1(1:length(ecgt),i+1)=ecgt;
% %     figure;
% %     plot(fftshift(abs(fft(ecgt1(:,i+1)))));
% end
% 
% for i=1:10
%     temp1=0;
%     temp2=0;
%     for mm=length(ecgt1(:,i)):-1:1
%         if(ecgt1(mm,i)~=0)
%             temp2=temp2+1;
%             temp1(temp2)=ecgt1(mm,i);
%         end
%     end
%     xyy=fliplr(temp1);
% %     plot(xyy);
%     afft=fftshift(abs(fft(xyy)));
%     figure(1)
%     subplot(2,5,i)
%     plot(xyy);
%     str=sprintf('%d. T',i*30+10);
%     title(str);
%     figure(2)
%     subplot(2,5,i)
%     plot(afft);
%     str=sprintf('%d. T için',i*30+10);
%     title(str);    
% end
% % Her bir fft i subplot etmek için deðerleri bir deðiþkene atmam lazým. fft
% % lerin uzunluðu 54 
% 
% % a=zeros(length(ecgt1),10);
% % 
% % for i=1:10
% %     a(:,i)=fftshift(abs(fft(ecgt1(:,i))));
% % end
% % FFT leri Çizdirmek için döngü
% % for i=0:9
% %     subplot(2,5,i+1)
% %     plot(a(:,i+1));
% %     str=sprintf('%d. T için',i*30+10);
% %     title(str);
% % end
% % 
% % % Sinyalleri çizdirmek için döngü
% % 
% % for i=1:10
% %     figure;
% %     plot(ecgt1(:,i));
% % end
% 
% else  
%     
% ecgt=0;
% ecgt1=0;
% % for i=1:8
% %     ecgt=ecgsignal(t_peak_start_pos(i):t_peak_final_pos(i));
% %     ecgt1(1:length(ecgt),i)=ecgt;
% %     figure;
% %     plot(fftshift(abs(fft(ecgt1(:,i)))));
% % end
% 
% % Burda ECG sinyalinin fft lerin içizdiriyorum. 10-40-70 etc.
% % sinyal çizdirme için
% 
% for i=0:7
%     ecgt=ecgfiltered_last(t_peak_start_pos(i*30+10):t_peak_final_pos(i*30+10));
%     ecgt1(1:length(ecgt),i+1)=ecgt;
% %     figure;
% %     plot(fftshift(abs(fft(ecgt1(:,i+1)))));
% end
% 
% for i=1:8
%     temp1=0;
%     temp2=0;
%     for mm=length(ecgt1(:,i)):-1:1
%         if(ecgt1(mm,i)~=0)
%             temp2=temp2+1;
%             temp1(temp2)=ecgt1(mm,i);
%         end
%     end
%     xyy=fliplr(temp1);
% %     plot(xyy);
%     afft=fftshift(abs(fft(xyy)));
%     figure(1)
%     subplot(2,5,i)
%     plot(xyy);
%     str=sprintf('%d. T',i*30+10);
%     title(str);
%     figure(2)
%     subplot(2,5,i)
%     plot(afft);
%     str=sprintf('%d. T için',i*30+10);
%     title(str);    
% end
% end

% Yeni kodlar

if(length(t_peak_start_pos)>250) 
ecgt=0;
ecgt1=0;
% for i=1:8
%     ecgt=ecgsignal(t_peak_start_pos(i):t_peak_final_pos(i));
%     ecgt1(1:length(ecgt),i)=ecgt;
%     figure;
%     plot(fftshift(abs(fft(ecgt1(:,i)))));
% end

% Burda ECG sinyalinin fft lerin içizdiriyorum. 10-40-70 etc.
% sinyal çizdirme için

Fx_eksen = linspace(-500,500-1000/512,512);  % Temel frekansý bulduktan sonra genlik olarak yarýsýna denk gelen nokta
                                               % bize bant geniþliðini
                                               % verir. Full width half
                                               % maximum bulunan frekansýn
                                               % 2 katý olarka yazýlacak
                                               % tabloya. lineer
                                               % interpolasyonla katsayýyý
                                               % bulup frekansa uygula,
                                               % frekans deðerini al. Sonra
                                               % 2 katýný al. Tabloya þu
                                               % þekilde yapcaz
                                               % frekans genliði (ex. 41.02
                                               % Hz) Bant geniþliði, FWHM de 2 katý olan deðer. 

for i=0:25
    ecgt=ecgfiltered_last(t_peak_start_pos(i*10+10):t_peak_final_pos(i*10+10));
    ecgt1(1:length(ecgt),i+1)=ecgt;
%     figure;
%     plot(fftshift(abs(fft(ecgt1(:,i+1)))));
end

band_freq_val=zeros(1,26);
fwhm=zeros(1,26);
frekansdegeri_xpoz=zeros(1,26);
max_genlik=zeros(1,26);
ortalamasure=zeros(1,26);

for i=1:26
    temp1=0;
    temp2=0;
    for mm=length(ecgt1(:,i)):-1:1
        if(ecgt1(mm,i)~=0)
            temp2=temp2+1;
            temp1(temp2)=ecgt1(mm,i);
        end
    end
    xyy=fliplr(temp1);
    xyyzeromean=xyy-mean(xyy);   % mean atma
%     plot(xyy);
    afft=fftshift(abs(fft(xyyzeromean,512)));  % FFT si sinyalin
    frekansdegeri=find(afft==max(afft(257:357))); % max frekansý bulma( temel frekans)
    frekansdegeri_max=frekansdegeri(2); % mirror olduðu için
    frekansdegeri_xeksen=Fx_eksen(frekansdegeri);  % frekans degerinin yerinin eksendeki karþýlýðý
    frekansdegeri_xpoz(i) = frekansdegeri_xeksen(2);   %frekans deðeri pozitif yeri
    max_genlik(i)=afft(frekansdegeri_max);     % Max frekans genliði
    band_genislik= max_genlik(i)/2;  % Bant genislik genligi
    
    val = band_genislik; %value to find
    degerler=afft(frekansdegeri_max:frekansdegeri_max+20);
    sonuclar=val<degerler;
    aaa=find(sonuclar==0);
    idx1=aaa(1);
%     [idx1 idx1] = min(temp3); %index of closest value
    cl1_eks=frekansdegeri_max+idx1-1;
    cl2_eks=cl1_eks-1;
    closest1 = afft(cl1_eks); %closest value 2
    closest2 = afft(cl2_eks); %closest value 2
    
    closest1_frekansval = Fx_eksen(cl1_eks);
    closest2_frekansval = Fx_eksen(cl2_eks);
    
    band_freq_val(i) = (((closest2_frekansval-closest1_frekansval)*(band_genislik-closest2)) / (closest2-closest1)) + closest2_frekansval;
    % frekans bant geniþliði
    bandpoz=find(Fx_eksen==closest1_frekansval);
    sondeger=find(afft==min(afft(bandpoz:bandpoz+12)));
    bbb=sondeger(2);
    fwhm(i)=Fx_eksen(bbb);
    ortalamasure(i)=length(xyyzeromean)/1000;

    
    figure(1)
    subplot(2,13,i)
    plot(xyyzeromean);
    str=sprintf('%d. T',i*10+10);
    title(str);
%     figure(2)
%     subplot(2,13,i)
%     plot(Fx_eksen,afft);
%     str=sprintf('%d. T için',i*10+10);
%     title(str);    
end
    basefilename='ECG_Grup5';
    dlmwrite([basefilename num2str(FileName([10 11 12 13])) '.csv'] ,[frekansdegeri_xpoz' max_genlik' band_freq_val' ortalamasure'],'precision','%20.5f');
    % Temel Frek(max), Bant Gen, Full width half max.


% Her bir fft i subplot etmek için deðerleri bir deðiþkene atmam lazým. fft
% lerin uzunluðu 54 

% a=zeros(length(ecgt1),10);
% 
% for i=1:10
%     a(:,i)=fftshift(abs(fft(ecgt1(:,i))));
% end
% FFT leri Çizdirmek için döngü
% for i=0:9
%     subplot(2,5,i+1)
%     plot(a(:,i+1));
%     str=sprintf('%d. T için',i*30+10);
%     title(str);
% end
% 
% % Sinyalleri çizdirmek için döngü
% 
% for i=1:10
%     figure;
%     plot(ecgt1(:,i));
% end

else  
    
ecgt=0;
ecgt1=0;

Fx_eksen = linspace(-500,500-1000/512,512);  % Temel frekansý bulduktan sonra genlik olarak yarýsýna denk gelen nokta

% for i=1:8
%     ecgt=ecgsignal(t_peak_start_pos(i):t_peak_final_pos(i));
%     ecgt1(1:length(ecgt),i)=ecgt;
%     figure;
%     plot(fftshift(abs(fft(ecgt1(:,i)))));
% end

% Burda ECG sinyalinin fft lerin içizdiriyorum. 10-40-70 etc.
% sinyal çizdirme için

for i=0:15
    ecgt=ecgfiltered_last(t_peak_start_pos(i*10+10):t_peak_final_pos(i*10+10));
    ecgt1(1:length(ecgt),i+1)=ecgt;
%     figure;
%     plot(fftshift(abs(fft(ecgt1(:,i+1)))));
end

band_freq_val=zeros(1,16);
fwhm=zeros(1,16);
frekansdegeri_xpoz=zeros(1,16);
max_genlik=zeros(1,16);
ortalamasure=zeros(1,16);

for i=1:16
    temp1=0;
    temp2=0;
    for mm=length(ecgt1(:,i)):-1:1
        if(ecgt1(mm,i)~=0)
            temp2=temp2+1;
            temp1(temp2)=ecgt1(mm,i);
        end
    end
%     xyy=fliplr(temp1);
%     xyyzeromean=xyy-mean(xyy);
% %     plot(xyy);
%     afft=fftshift(abs(fft(xyyzeromean,32)));   %%%%% FFT boyutunu hepsini ayný alalým.
%     figure(1)
%     subplot(2,5,i)
%     plot(xyy-mean(xyy));
%     str=sprintf('%d. T',i*30+10);
%     title(str);
%     figure(2)
%     subplot(2,5,i)
%     plot(afft);
%     str=sprintf('%d. T için',i*30+10);
%     title(str);    
    xyy=fliplr(temp1);
    xyyzeromean=xyy-mean(xyy);   % mean atma
%     plot(xyy);
    afft=fftshift(abs(fft(xyyzeromean,512)));  % FFT si sinyalin
    frekansdegeri=find(afft==max(afft(257:357))); % max frekansý bulma( temel frekans)
    frekansdegeri_max=frekansdegeri(2); % mirror olduðu için
    frekansdegeri_xeksen=Fx_eksen(frekansdegeri);  % frekans degerinin yerinin eksendeki karþýlýðý
    frekansdegeri_xpoz(i) = frekansdegeri_xeksen(2);   %frekans deðeri pozitif yeri
    max_genlik(i)=afft(frekansdegeri_max);     % Max frekans genliði
    band_genislik= max_genlik(i)/2;  % Bant genislik genligi
    
    val = band_genislik; %value to find
    degerler=afft(frekansdegeri_max:frekansdegeri_max+20);
    sonuclar=val<degerler;
    aaa=find(sonuclar==0);
    idx1=aaa(1);
%     [idx1 idx1] = min(temp3); %index of closest value
    cl1_eks=frekansdegeri_max+idx1-1;
    cl2_eks=cl1_eks-1;
    closest1 = afft(cl1_eks); %closest value 2
    closest2 = afft(cl2_eks); %closest value 2
    
    closest1_frekansval = Fx_eksen(cl1_eks);
    closest2_frekansval = Fx_eksen(cl2_eks);
    
    band_freq_val(i) = (((closest2_frekansval-closest1_frekansval)*(band_genislik-closest2)) / (closest2-closest1)) + closest2_frekansval;
    % frekans bant geniþliði
    bandpoz=find(Fx_eksen==closest1_frekansval);
    sondeger=find(afft==min(afft(bandpoz:bandpoz+12)));
    bbb=sondeger(2);
    fwhm(i)=Fx_eksen(bbb);
    ortalamasure(i)=length(xyyzeromean)/1000;

    basefilename='ECG_Grup5';
    dlmwrite([basefilename num2str(FileName([10 11 12 13])) '.csv'] ,[frekansdegeri_xpoz' max_genlik' band_freq_val' ortalamasure'],'precision','%20.5f');
    % Temel Frek(max), Bant Gen, Full width half max.

    figure(1)
    subplot(2,8,i)
    plot(xyyzeromean);
    str=sprintf('%d. T',i*10+10);
    title(str);
%     figure(2)
%     subplot(2,5,i)
%     plot(Fx_eksen,afft);
%     str=sprintf('%d. T için',i*30+10);
%     title(str);    
end
end

figure
plot(ecgfiltered_last)