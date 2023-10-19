clear all;
close all;

subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

path=['D:\jiangxinyu\Open Access HDsEMG DataSet Analysis\hyser_dataset']; % change path to the location of the dataset

zc_ssc_thresh=0.0004; %threshold to detect valid zero cross and slope sign change features
window_len=0.02; % 0.02 s sliding window to extract features
step_len=0.02; % 0.02 s step length of the sliding window to extract features
f_cutoff=10; % filter force data using a low pass filter with 10 Hz cut off frequency
fs_emg=2048;
fs_force=100;
Q=20;
D=1;
ii=0;
Tol=0.05;
pca_active=1;
dim=200;

parfor i=1:20
    for session=1:2
        
        i
        session
        
        force_finger_random=load_random(path,subject{1,i},session,'force');
        emg_finger_random=load_random(path,subject{1,i},session,'preprocess');
        
        mvc=get_mvc(path,subject{1,i},num2str(session));
        force_norm=normalize_force(force_finger_random,mvc);
        force_norm_preprocess=preprocess_force(force_norm,window_len,step_len,f_cutoff,fs_force,fs_emg);
        
        feature=cell(1,5);
        force_groundtruth_thumb=cell(1,5);
        force_groundtruth_index=cell(1,5);
        force_groundtruth_middle=cell(1,5);
        force_groundtruth_ring=cell(1,5);
        force_groundtruth_little=cell(1,5);
        
        for u=1:size(emg_finger_random,2)
            emg=emg_finger_random{1,u};
            rms=get_rms(emg,window_len,step_len,fs_emg);
            wl=get_wl(emg,window_len,step_len,fs_emg);
            zc=get_zc(emg,window_len,step_len,zc_ssc_thresh, fs_emg);
            ssc=get_ssc(emg,window_len,step_len,zc_ssc_thresh, fs_emg);
            feature{1,u}=[rms';wl';zc';ssc'];
            force_groundtruth_thumb{1,u}=force_norm_preprocess{1,u}(:,1);
            force_groundtruth_index{1,u}=force_norm_preprocess{1,u}(:,2);
            force_groundtruth_middle{1,u}=force_norm_preprocess{1,u}(:,3);
            force_groundtruth_ring{1,u}=force_norm_preprocess{1,u}(:,4);
            force_groundtruth_little{1,u}=force_norm_preprocess{1,u}(:,5);
        end
        
        % model training and testing
        for test_trial_idx=1:5
            
            force_test_thumb=force_groundtruth_thumb(1,test_trial_idx);
            force_train_thumb=force_groundtruth_thumb;
            force_train_thumb(:,test_trial_idx)=[];
            
            force_test_index=force_groundtruth_index(1,test_trial_idx);
            force_train_index=force_groundtruth_index;
            force_train_index(:,test_trial_idx)=[];
            
            force_test_middle=force_groundtruth_middle(1,test_trial_idx);
            force_train_middle=force_groundtruth_middle;
            force_train_middle(:,test_trial_idx)=[];
            
            force_test_ring=force_groundtruth_ring(1,test_trial_idx);
            force_train_ring=force_groundtruth_ring;
            force_train_ring(:,test_trial_idx)=[];
            
            force_test_little=force_groundtruth_little(1,test_trial_idx);
            force_train_little=force_groundtruth_little;
            force_train_little(:,test_trial_idx)=[];
                
            feature_test=feature(1,test_trial_idx);
            feature_train=feature;
            feature_train(:,test_trial_idx)=[];
   
            [feature_train_norm,feature_test_norm]=feature_normalize(feature_train,feature_test,pca_active,dim);
                
            % train the regression model 
            x_plus_thumb = e_sifir_trn(feature_train_norm, force_train_thumb, [Q D Tol ii], [2/window_len 2/window_len]);
            x_plus_index = e_sifir_trn(feature_train_norm, force_train_index, [Q D Tol ii], [2/window_len 2/window_len]);
            x_plus_middle = e_sifir_trn(feature_train_norm, force_train_middle, [Q D Tol ii], [2/window_len 2/window_len]);
            x_plus_ring = e_sifir_trn(feature_train_norm, force_train_ring, [Q D Tol ii], [2/window_len 2/window_len]);
            x_plus_little = e_sifir_trn(feature_train_norm, force_train_little, [Q D Tol ii], [2/window_len 2/window_len]);
            
            [Err_thumb(test_trial_idx,i,session), force_estimate_thumb] = e_sifir_tst(x_plus_thumb, feature_test_norm, force_test_thumb, [Q D Tol ii], [2/window_len 2/window_len]);           
            [Err_index(test_trial_idx,i,session), force_estimate_index] = e_sifir_tst(x_plus_index, feature_test_norm, force_test_index, [Q D Tol ii], [2/window_len 2/window_len]);               
            [Err_middle(test_trial_idx,i,session), force_estimate_middle] = e_sifir_tst(x_plus_middle, feature_test_norm, force_test_middle, [Q D Tol ii], [2/window_len 2/window_len]);               
            [Err_ring(test_trial_idx,i,session), force_estimate_ring] = e_sifir_tst(x_plus_ring, feature_test_norm, force_test_ring, [Q D Tol ii], [2/window_len 2/window_len]);               
            [Err_little(test_trial_idx,i,session), force_estimate_little] = e_sifir_tst(x_plus_little, feature_test_norm, force_test_little, [Q D Tol ii], [2/window_len 2/window_len]);               
            
        end
    end
end
average_error_thumb=mean(mean(mean(Err_thumb,3),2),1)% unit of error: %MVC 
average_error_index=mean(mean(mean(Err_index,3),2),1)
average_error_middle=mean(mean(mean(Err_middle,3),2),1)
average_error_ring=mean(mean(mean(Err_ring,3),2),1)
average_error_little=mean(mean(mean(Err_little,3),2),1)

save('../results/Err_thumb_random.mat','Err_thumb');
save('../results/Err_index_random.mat','Err_index');
save('../results/Err_middle_random.mat','Err_middle');
save('../results/Err_ring_random.mat','Err_ring');
save('../results/Err_little_random.mat','Err_little');