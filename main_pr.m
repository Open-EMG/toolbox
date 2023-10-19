clear all;
close all;

subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

path=['D:\jiangxinyu\Open Access HDsEMG DataSet Analysis\hyser_dataset']; % change path to the location of the dataset

task_type='dynamic'; % change to 'maintenance' if you need to classify maintenance tasks

switch task_type
    case 'dynamic'
        window_len=0.75;
        step_len=0.75;
    case 'maintenance'
        window_len=3.75;
        step_len=3.75;
end
sig_start=0.25;% remove the first 0.25s startup duration.

zc_ssc_thresh=0.0004; %threshold to detect valid zero cross and slope sign change features
fs_emg=2048;
pca_active=1;

parfor i=1:20
    for session=1:2
        
        i
        session
        
        [data,label]=load_pr(path,subject{1,i},session,task_type,'preprocess');
        
        Ntrial=length(data);
        feature=[];
        for j=1:Ntrial
            emg=data{1,j}(sig_start*fs_emg+1:end,:);
            rms_tmp=get_rms(emg,window_len,step_len,fs_emg);
            rms=reshape(rms_tmp,[1,numel(rms_tmp)]);
            wl_tmp=get_wl(emg,window_len,step_len,fs_emg);
            wl=reshape(wl_tmp,[1,numel(wl_tmp)]);
            zc_tmp=get_zc(emg,window_len,step_len,zc_ssc_thresh,fs_emg);
            zc=reshape(zc_tmp,[1,numel(zc_tmp)]);
            ssc_tmp=get_ssc(emg,window_len,step_len,zc_ssc_thresh,fs_emg);
            ssc=reshape(ssc_tmp,[1,numel(ssc_tmp)]);
            feature(:,j)=[rms';wl';zc';ssc'];
        end
        
        predict_label=zeros(1,Ntrial);
        for j=1:Ntrial
            feature_test=feature(:,j);
            feature_train=feature;
            feature_train(:,j)=[];
            
            label_test=label(j);
            label_train=label;
            label_train(j)=[];
            
            dim=size(feature_train,2)-1;
            [feature_train_norm,feature_test_norm]=feature_normalize(feature_train,feature_test,pca_active,dim);
            
            mdl = ClassificationDiscriminant.fit(feature_train_norm',label_train);
            predict_label(j) = predict(mdl, feature_test_norm');
        end
        accuracy(i,session)=mean(double((predict_label==label)));
    end
end
average_accuracy=mean(mean(accuracy))  
save(['../results/accuracy_',task_type,'.mat'],'accuracy');