function [EMGOutput,W] = SimEMGProcessing(EMGInput,varargin)
% Set Default
SNR = 20;
WhitenFlag = 'On';
R = 16;
NeedleType = 'MonoPolar';

for i = 1:2:length(varargin)
    switch varargin{i}
        
        case 'SNR'
            SNR = varargin{i+1};
        case 'WhitenFlag'
            WhitenFlag = varargin{i+1};
        case 'R'
            R = varargin{i+1};
        case 'NeedleType'
            NeedleType = varargin{i+1};
    end
end
[N,NumCh1] = size(EMGInput);
if strcmp(SNR,'Inf')==1
    EMGNoisedtemp = EMGInput;
else
    EMGNoisedtemp = awgn(EMGInput,SNR,'measured');
end

if strcmp(NeedleType,'MonoPolar')==1
    EMGNoised  = EMGNoisedtemp ;
else if strcmp(NeedleType,'BiPolar')==1
        EMGNoised = zeros(N,NumCh1-8);
        for i = 1:NumCh1-8
            EMGNoised(:,i) = EMGNoisedtemp(:,i) - EMGNoisedtemp(:,i+8);
        end
    end
end
[~,NumCh2] = size(EMGNoised);
EMGExtended = zeros(N,NumCh2*(R+1));
EMGExtended(:,1:NumCh2) = EMGNoised;
if R~=0
    for i = 1:R
        EMGExtended(1+i:end,NumCh2*i+1:NumCh2*i+NumCh2) = EMGNoised(1:end-i,:);
    end
end
EMGExtendedSubMean = EMGExtended - ones(N,1)*mean(EMGExtended,1);

if strcmp(WhitenFlag,'On')==1
    Rxx = EMGExtendedSubMean'*EMGExtendedSubMean;
    [V,D] = eig(Rxx);
    fudgefactor = 0;
%     W = (V*diag(1./(diag(D)+fudgefactor).^(1/2))*V')';
    W = real((V*diag(1./(diag(D)+fudgefactor).^(1/2))*V'))';
    EMGOutput = W * EMGExtendedSubMean';
else
    EMGOutput = EMGExtendedSubMean';
    W = [];
end
    
end

