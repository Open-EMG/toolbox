clear all;
close all;

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

subject='01';
session=1;
finger=2;%index finger

force_finger_1dof=load_1dof(path,subject,session,'force');
emg_finger_1dof=load_1dof(path,subject,session,'preprocess');

% take all 3 samples of finger 2 acquired in 3 trials;
force_finger_1dof=force_finger_1dof(finger,:);
emg_finger_1dof=emg_finger_1dof(finger,:);
        
mvc=get_mvc(path,subject,num2str(session));
force_norm=normalize_force(force_finger_1dof,mvc);
force_norm_preprocess=preprocess_force(force_norm,window_len,step_len,f_cutoff,fs_force,fs_emg);

feature=cell(1,3);
force_groundtruth=cell(1,3);
for v=1:3
    emg=emg_finger_1dof{1,v};
    rms=get_rms(emg,window_len,step_len,fs_emg);
    wl=get_wl(emg,window_len,step_len,fs_emg);
    zc=get_zc(emg,window_len,step_len,thresh, fs_emg);
    ssc=get_ssc(emg,window_len,step_len,thresh, fs_emg);
    feature{1,v}=[rms';wl';zc';ssc'];
    force_groundtruth{1,v}=force_norm_preprocess{1,v}(:,finger);
end

% model training and testing     
train_trial_idx=1 % use the 1st trial as the training data
test_trial_idx=2 % use the 2nd trial as the testing data
force_train=force_groundtruth(1,train_trial_idx);
force_test=force_groundtruth(1,test_trial_idx);
  
feature_train=feature(1,train_trial_idx);
feature_test=feature(1,test_trial_idx);
   
[feature_train_norm,feature_test_norm]=feature_normalize(feature_train,feature_test,pca_active,dim);
            
% train the regression model 
x_plus = e_sifir_trn(feature_train_norm, force_train, [Q D Tol ii], [2/window_len 2/window_len]);

[Err,force_estimate] = e_sifir_tst(x_plus, feature_test_norm, force_test, [Q D Tol ii], [2/window_len 2/window_len]);
