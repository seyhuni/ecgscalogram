for i=1:5
[info,data]=acqread( ['Sevgi Karaboyun-4.acq']);
ECG=double(data{1,2});
eval(['Sevgi Karaboyun-4' num2str(i) '=' 'ECG' ';']); 
save(['Sevgi Karaboyun-4' num2str(i)],['Sevgi Karaboyun-4' num2str(i)]);
clear data info ECG
end

[info,data]=acqread('TEB-MCT+BOS-A5-y.acq');
ECG=double(data{1,2});
ECG_TEB_MCT_BOS_A5_y = ECG; 
save('ECG-TEB-MCT+BOS-A5-y' , 'ECG_TEB_MCT_BOS_A5_y');
clear data info ECG

for i=1:5
[info,data]=acqread( ['TEB-MCT+BOS-B' num2str(i) '.acq']);
ECG=double(data{1,2});
eval(['ECG_TEB_MCT_BOS_B' num2str(i) '=' 'ECG' ';']); 
save(['ECG-TEB-MCT+BOS-B' num2str(i)],['ECG_TEB_MCT_BOS_B' num2str(i)]);
clear data info ECG
end

[info,data]=acqread('TEB-MCT+BOS-B2-y.acq');
ECG=double(data{1,2});
ECG_TEB_MCT_BOS_B2_y = ECG; 
save('ECG-TEB-MCT+BOS-B2-y' , 'ECG_TEB_MCT_BOS_B2_y');
clear data info ECG

[info,data]=acqread('TEB-MCT+BOS-B2-y1.acq');
ECG=double(data{1,2});
ECG_TEB_MCT_BOS_B2_y1 = ECG; 
save('ECG-TEB-MCT+BOS-B2-y1' , 'ECG_TEB_MCT_BOS_B2_y1');
clear data info ECG

[info,data]=acqread('TEB-MCT+BOS-B5-y1.acq');
ECG=double(data{1,2});
ECG_TEB_MCT_BOS_B5_y1 = ECG; 
save('ECG-TEB-MCT+BOS-B5-y1' , 'ECG_TEB_MCT_BOS_B5_y1');
clear data info ECG





% Ts=1e-3;
% data=ECG;
% t=0:Ts:Ts*(length(data)-1);
% figure; plot(t,data); grid on
% % data1=data(1400/Ts:4600/Ts);
% % plot(linspace(1400,4600,length(data1)),data1); grid on
% % data1=data1-mean(data1);
% f=(0:length(data)-1)./length(data);
% f=f.*1000;
% fftdata=abs(fft(data));
% fftdata=2*fftdata/length(fftdata);
% plot(f(1:length(f)/2+1),fftdata(1:length(f)/2+1)); grid on
% semilogx(f(1:length(f)/2+1),fftdata(1:length(f)/2+1)); grid on