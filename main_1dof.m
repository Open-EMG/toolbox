clear all;
close all;

subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

path=['D:\jiangxinyu\Open Access HDsEMG DataSet Analysis\hyser_dataset']; % change path to the location of the dataset

thresh=0.0004 % 0.4 mv threshold to detect valid zero cross and slope sign change features
window_len=0.02 % 0.02 s sliding window to extract features
step_len=0.02 % 0.02 s step length of the sliding window to extract features
f_cutoff=10; % filter force data using a low pass filter with 10 Hz cut off frequency
fs_emg=2048;
fs_force=100;
Q=20;
ii=0;
D=1;
Tol=0.05;
pca_active=1;
dim=200;

parfor i=1:20
    for session=1:2
        
        i
        session
        
        force_finger_1dof=load_1dof(path,subject{1,i},session,'force');
        emg_finger_1dof=load_1dof(path,subject{1,i},session,'preprocess');
        
        mvc=get_mvc(path,subject{1,i},num2str(session));
        force_norm=normalize_force(force_finger_1dof,mvc);
        force_norm_preprocess=preprocess_force(force_norm,window_len,step_len,f_cutoff,fs_force,fs_emg);

        feature=cell(5,3);
        force_groundtruth=cell(5,3);
        for u=1:5
            for v=1:3
                emg=emg_finger_1dof{u,v};
                rms=get_rms(emg,window_len,step_len,fs_emg);
                wl=get_wl(emg,window_len,step_len,fs_emg);
                zc=get_zc(emg,window_len,step_len,thresh, fs_emg);
                ssc=get_ssc(emg,window_len,step_len,thresh, fs_emg);
                feature{u,v}=[rms';wl';zc';ssc'];
                force_groundtruth{u,v}=force_norm_preprocess{u,v}(:,u);
            end
        end
        
        % model training and testing
        for test_finger_idx=1:5
            for test_trial_idx=1:3
                
                force_test=force_groundtruth(test_finger_idx,test_trial_idx);
                force_train=force_groundtruth(test_finger_idx,:);
                force_train(:,test_trial_idx)=[];
                
                feature_test=feature(test_finger_idx,test_trial_idx);
                feature_train=feature(test_finger_idx,:);
                feature_train(:,test_trial_idx)=[];
   
                [feature_train_norm,feature_test_norm]=feature_normalize(feature_train,feature_test,pca_active,dim);
            
                % train the regression model 
                x_plus = e_sifir_trn(feature_train_norm, force_train, [Q D Tol ii], [2/window_len 2/window_len]);

                [Err(test_finger_idx,test_trial_idx,i,session),force_estimate] = e_sifir_tst(x_plus, feature_test_norm, force_test, [Q D Tol ii], [2/window_len 2/window_len]);
                
            end
        end
    end
end
mean(mean(mean(mean(Err,4),3),2),1)  
save('../results/Err_1dof.mat','Err');