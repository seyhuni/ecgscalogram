srcFiles = dir('C:\Users\huseyinyanik\Desktop\akademik\EKG ÇALIÞMALARI\antrene work\kayitlar\Kontrol Grubu Kayýtlarý\21.04.2011-yeni kontrol\Erkek Kontrol\*.acq');
for i=1:length(srcFiles)
filename = strcat('C:\Users\huseyinyanik\Desktop\akademik\EKG ÇALIÞMALARI\antrene work\kayitlar\Kontrol Grubu Kayýtlarý\21.04.2011-yeni kontrol\Erkek Kontrol\',srcFiles(i).name)
[info,data]=acqread(filename);
ECG=double(data{1,2});
save([srcFiles(i).name(1:end-4)]);
clear data info ECG
end