function [FiringRate,TimeWin] = FiringRateEst(SpikeTrain,GroupIndex,WinLength,StepLength,FsEMG,MeanFlag)
[NumGroup,~] = size(GroupIndex);
[N,~] = size(SpikeTrain);

WinNumSample = floor(FsEMG*WinLength);
StepNumSampe = floor(FsEMG*StepLength);
count = 1;
for j = 1:StepNumSampe:N-WinNumSample+1
    for ii = 1:NumGroup
        Indextemp = nonzeros(GroupIndex(ii,:));
        SpikeSegtemp = SpikeTrain(j:WinNumSample+j-1,Indextemp);
        Firingtemp = sum(SpikeSegtemp,1);
        FiringNumMU1temp = length(find(Firingtemp~=0));
        if FiringNumMU1temp~=0
            FiringRate(count,ii) = sum(sum(SpikeSegtemp))/WinLength;
        else
            FiringRate(count,ii) = 0;
        end
    end    
    count = count+1;
end
TimeMax = (j+WinNumSample/2)/FsEMG;
TimeWin = WinNumSample/2/FsEMG:StepNumSampe/FsEMG:TimeMax;

if strcmp(MeanFlag,'On')
    for ii = 1:NumGroup
        Indextemp = nonzeros(GroupIndex(ii,:));
        FiringRate = FiringRate/length(Indextemp);
    end
end
end