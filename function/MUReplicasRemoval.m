function [SpikeTrain,Goodindex] = MUReplicasRemoval(SpikeTrain,s1,Fs)
% Post-process the results of sEMG decomposition via physiological basis. 
% Goodindex returns the indices of motor units which are not noises or
% motion artifacts.

Timetemp = (1/Fs:1/Fs:length(SpikeTrain)/Fs)';

% Step 1
Firings = sum(SpikeTrain,1);
index1 = find(Firings>4*Timetemp(end));
index2 = find(Firings<35*Timetemp(end));
Goodindextemp = intersect(index1,index2);
NumGood = length(Goodindextemp);

% Step 2
Time = Timetemp*ones(1,NumGood);
FirT = cell(NumGood,1);
for k = 1:NumGood
    loc = find(SpikeTrain(:,Goodindextemp(k))==1);
    Diffloc = diff(loc);
    loc2 = Diffloc<Fs*0.02;
    for l = 1:length(loc2)
        if loc2(l) == 1
           peaktemp1 = s1(loc(l),k);
           peaktemp2 = s1(loc(l+1),k);
           if peaktemp1>=peaktemp2
           SpikeTrain(loc(l+1),Goodindextemp(k)) = 0;
           else
               SpikeTrain(loc(l),Goodindextemp(k)) = 0;
           end
        end
    end
end

for k = 1:NumGood
    loc = find(SpikeTrain(:,Goodindextemp(k))==1);
    Diffloc = diff(loc);
    loc2 = Diffloc<Fs*0.02;
    for l = 1:length(loc2)
        if loc2(l) == 1
           peaktemp1 = s1(loc(l),k);
           peaktemp2 = s1(loc(l+1),k);
           if peaktemp1>=peaktemp2
           SpikeTrain(loc(l+1),Goodindextemp(k)) = 0;
           else
               SpikeTrain(loc(l),Goodindextemp(k)) = 0;
           end
        end
    end
end

for k = 1:NumGood
    loc = find(SpikeTrain(:,Goodindextemp(k))==1);
    Diffloc = diff(loc);
    loc2 = Diffloc<Fs*0.02;
    for l = 1:length(loc2)
        if loc2(l) == 1
           peaktemp1 = s1(loc(l),k);
           peaktemp2 = s1(loc(l+1),k);
           if peaktemp1>=peaktemp2
           SpikeTrain(loc(l+1),Goodindextemp(k)) = 0;
           else
               SpikeTrain(loc(l),Goodindextemp(k)) = 0;
           end
        end
    end
    FirT{k} = Time(SpikeTrain(:,Goodindextemp(k))==1);
end

% Step 3
NumMU = length(FirT);
count = 1;
index = 1:NumMU;
NumRemoval = 0;

wrong=0;
while length(index)~=count
    indexRemovaltemp = [];
    for i = 1:length(index)-count-1
        Logic = CSIndex(FirT{count}, FirT{(count+i)}, 0.01, 10);
        if Logic == 1
            indexRemovaltemp = [indexRemovaltemp count+i];
        end
    end
    FirT(indexRemovaltemp) = [];
    index(indexRemovaltemp) =[];
    NumRemoval = length(indexRemovaltemp);
    count = count+1;
     if(count>10000)
        wrong=1;
        break;
    end
end
Goodindex= Goodindextemp(index);

if(wrong==1)
    Goodindex=[];
    SpikeTrain=[];
end