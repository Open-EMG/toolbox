function SIL = SILCal(s,Fs)
[b,a] = butter(4,500/(Fs/2),'low');
[~,NumMU] = size(s);
SIL = zeros(1,NumMU);
for i =1:NumMU

s(:,i) = filtfilt(b,a,s(:,i));
[pks,~] = findpeaks(s(:,i).^2);

[idx,~] = kmeansplus(pks',2);
% [idx,~] = myCluster2(pks);

sil = silhouette(pks,idx);
SIL(i) = (mean(sil(idx==1))+mean(sil(idx==2)))/2;
end
