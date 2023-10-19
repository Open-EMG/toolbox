function [s,B,SpikeTrain] = FastICA(EMG,M,Fs)
% EMG (N * M): M is the length of signal. N is the number of channels
% M: number of total iteration
% Fs: sample frequency
% s: separated source after decomposition
% B: separation vector
% SpikeTrain: motor unit spike trains after decomposition
% Default Settings
Tolx = 10^-4;
[NumCh,N] = size(EMG);
s = zeros(N,M);
stemp = zeros(N,M);
B = zeros(NumCh,1);
SpikeTrain = zeros(N,M);

% W = eye(NumCh);
% [b,a] = butter(4,100/(Fs/2),'low');
% Gx = x^3/3; gx = x^2;; gx' = 2x;
for i = 1:M
    w = [];
    %     w(:,1) = zeros(NumCh,1);
    %     w(:,2) = W(:,i);
    w(:,1) = randn(NumCh,1);
    w(:,2) = randn(NumCh,1);
    for n = 2:100
        if abs(w(:,n)'*w(:,n-1)-1)>Tolx
            A = mean(2*w(:,n)'*EMG);
            w(:,n+1) = EMG*(((w(:,n)'*EMG)').^2)-A*w(:,n);
            w(:,n+1) = w(:,n+1) - B*B'*w(:,n+1);
            w(:,n+1) = w(:,n+1)/norm(w(:,n+1));
        else
            break;
        end
    end
    CoV(1) = 1;
    CoV(2) = 0.99;
    %     for m = n:100
    %         if abs(CoV(m-n+2)-CoV(m-n+1))>Tolx
%     stemp(:,i) = w(:,n)'*EMG;
%     s(:,i) = filtfilt(b,a,stemp(:,i));
    s(:,i) = w(:,n)'*EMG;
    [pks,loc] = findpeaks(s(:,i).^2);
    
    
    [idx,C] = kmeansplus(pks',2);
%     [idx,~] = myCluster2(pks);
    
    
    if sum(idx==1)<=sum(idx==2)
        SpikeLoc = loc(idx==1);
    else
        SpikeLoc = loc(idx==2);
    end
    %             SpikeInterval = diff(SpikeLoc);
    %             CoV(m-n+3) = std(SpikeInterval)/mean(SpikeInterval);
    %             w(:,m+1) = mean(EMG(:,SpikeLoc),2);
    %         else
    %             break;
    %         end
    %     end
    SpikeTrain(SpikeLoc,i) = 1;
    B(:,i) = w(:,end);
end