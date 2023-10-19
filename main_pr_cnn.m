clear all;
close all;

%% parameter setting
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

path=['D:\jiangxinyu\Open Access HDsEMG DataSet Analysis\hyser_dataset']; % change path to the location of the dataset

task_type='dynamic'; % change to 'maintenance' if you need to classify maintenance tasks

switch task_type
    case 'dynamic'
        K=6; %6-fold cross-validation
    case 'maintenance'
        K=2; %2-fold cross-validation
end
sig_start=0.25;% remove the first 0.25s startup duration.
step_len=0.125;
window_len=0.25;
fs_emg=2048;

% layout: convert 256 electrodes to 16 * 16 array
for i=1:8
    for j=1:8
        layout_array(i,j)=64-(i-1)*8-(j-1);
    end
end
layout=[[layout_array;64+layout_array],[128+layout_array;192+layout_array]];

%% data preparation

data_input=cell(2,K); % 2 sessions, K folds cross-validation in each session
label_input=cell(2,K);

for session=1:2 % index of experiment sessions
    for i=1:20 % index of subjects
        
        i
        session
        
        [data,label]=load_pr(path,subject{1,i},session,task_type,'preprocess');
        
        INDICES = crossvalind('Kfold',label,K); 
        
        win=hamming(256);
        for j=1:length(data)
            for t=floor(sig_start*fs_emg):floor(step_len*fs_emg):length(data{1,j})-floor(window_len*fs_emg)
                sig=data{1,j}(t+1:t+floor(window_len*fs_emg),:);
                [sig_stft,f,t]=stft(sig,fs_emg,'Window',win,'OverlapLength',128,'FrequencyRange','onesided');                               
                sig_tensor=reshape(abs(sig_stft(1:64,:,reshape(layout,[1,256]))),[64*3,16,16]); 
                sig_tensor_permute=permute(sig_tensor,[2 3 1]);
                if(~length(data_input{session,INDICES(j)}))
                    data_input{session,INDICES(j)}=sig_tensor_permute;
                    label_input{session,INDICES(j)}=label(j);
                else
                    data_input{session,INDICES(j)}(:,:,:,end+1)=sig_tensor_permute;
                    label_input{session,INDICES(j)}(end+1)=label(j);
                end
            end
        end
    end
end
    
%% define the network architecture

layers = [imageInputLayer([16,16,192])
          convolution2dLayer(3,32,'Padding','same')
          batchNormalizationLayer
          leakyReluLayer(0.1)
          dropoutLayer(0.2)
          maxPooling2dLayer(3,'Stride',2)
          convolution2dLayer(3,64,'Padding','same')
          batchNormalizationLayer
          leakyReluLayer(0.1)
          dropoutLayer(0.2)
          maxPooling2dLayer(3,'Stride',2)
          fullyConnectedLayer(576)
          fullyConnectedLayer(34)
          softmaxLayer
          classificationLayer];
      
options = trainingOptions('adam', ...
    'InitialLearnRate',0.001, ...
    'GradientDecayFactor',0.9, ...
    'SquaredGradientDecayFactor',0.999, ...
    'L2Regularization',0.0001, ...
    'MaxEpochs',50, ...
    'MiniBatchSize',512, ...
    'Plots','none')
 
for session=1:2
    for testFoldIdx=1:K
        data_test=data_input{session,testFoldIdx};
        label_test=label_input{session,testFoldIdx};
        train_fold=[1:K];
        train_fold(testFoldIdx)=[];
        data_train=[];
        label_train=[];
        for i=1:length(train_fold)
            data_train=cat(4,data_train,data_input{session,train_fold(i)});
            label_train=[label_train,label_input{session,train_fold(i)}];
        end
        convnet = trainNetwork(data_train,categorical(label_train),layers,options);
        [label_predict,score] = classify(convnet,data_test);
        
        Naverage=(size(data{1,1},1)-ceil(sig_start*fs_emg)-ceil(window_len*fs_emg))/ceil(step_len*fs_emg)+1;
        for i=1:length(label_test)/Naverage
            score_window_average(i,:)=mean(score((i-1)*Naverage+1:i*Naverage,:));
            label_predict_average(i)=find(score_window_average(i,:)==max(score_window_average(i,:)));
            label_test_average(i)=label_test(i*Naverage);
        end
        accuracy(session,testFoldIdx) = sum(label_predict_average == label_test_average)./numel(label_test_average)
    end
end
save(['../results/accuracy_',task_type,'_cnn.mat'],'accuracy');