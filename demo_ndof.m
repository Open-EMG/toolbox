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
combination=2;% finger combination (2)

force_finger_ndof=load_ndof(path,subject,session,'force');
emg_finger_ndof=load_ndof(path,subject,session,'preprocess');

force_finger_ndof=force_finger_ndof(combination,:);
emg_finger_ndof=emg_finger_ndof(combination,:);
        
mvc=get_mvc(path,subject,num2str(session));
force_norm=normalize_force(force_finger_ndof,mvc);
force_norm_preprocess=preprocess_force(force_norm,window_len,step_len,f_cutoff,fs_force,fs_emg);

feature=cell(1,2);
force_groundtruth_thumb=cell(1,2);
force_groundtruth_index=cell(1,2);
force_groundtruth_middle=cell(1,2);
force_groundtruth_ring=cell(1,2);
force_groundtruth_little=cell(1,2);
        
for v=1:size(emg_finger_ndof,2)
    v
    emg=emg_finger_ndof{1,v};
    rms=get_rms(emg,window_len,step_len,fs_emg);
    wl=get_wl(emg,window_len,step_len,fs_emg);
    zc=get_zc(emg,window_len,step_len,thresh, fs_emg);
    ssc=get_ssc(emg,window_len,step_len,thresh, fs_emg);
    feature{1,v}=[rms';wl';zc';ssc'];
    force_groundtruth_thumb{1,v}=force_norm_preprocess{1,v}(:,1);
    force_groundtruth_index{1,v}=force_norm_preprocess{1,v}(:,2);
    force_groundtruth_middle{1,v}=force_norm_preprocess{1,v}(:,3);
    force_groundtruth_ring{1,v}=force_norm_preprocess{1,v}(:,4);
    force_groundtruth_little{1,v}=force_norm_preprocess{1,v}(:,5);
end
        
% model training and testing     
train_trial_idx=1; % use the 1st trial as the training data
test_trial_idx=2; % use the 2nd trial as the testing data
force_train_thumb=force_groundtruth_thumb(1,train_trial_idx);
force_test_thumb=force_groundtruth_thumb(1,test_trial_idx);
force_train_index=force_groundtruth_index(1,train_trial_idx);
force_test_index=force_groundtruth_index(1,test_trial_idx);
force_train_middle=force_groundtruth_middle(1,train_trial_idx);
force_test_middle=force_groundtruth_middle(1,test_trial_idx);
force_train_ring=force_groundtruth_ring(1,train_trial_idx);
force_test_ring=force_groundtruth_ring(1,test_trial_idx);
force_train_little=force_groundtruth_little(1,train_trial_idx);
force_test_little=force_groundtruth_little(1,test_trial_idx);

feature_train=feature(1,train_trial_idx);
feature_test=feature(1,test_trial_idx);

[feature_train_norm,feature_test_norm]=feature_normalize(feature_train,feature_test,pca_active,dim);
            
% train the regression model 
x_plus_thumb = e_sifir_trn(feature_train_norm, force_train_thumb, [Q D Tol ii], [2/window_len 2/window_len]);
x_plus_index = e_sifir_trn(feature_train_norm, force_train_index, [Q D Tol ii], [2/window_len 2/window_len]);
x_plus_middle = e_sifir_trn(feature_train_norm, force_train_middle, [Q D Tol ii], [2/window_len 2/window_len]);
x_plus_ring = e_sifir_trn(feature_train_norm, force_train_ring, [Q D Tol ii], [2/window_len 2/window_len]);
x_plus_little = e_sifir_trn(feature_train_norm, force_train_little, [Q D Tol ii], [2/window_len 2/window_len]);
            
[Err_thumb, force_estimate_thumb] = e_sifir_tst(x_plus_thumb, feature_test_norm, force_test_thumb, [Q D Tol ii], [2/window_len 2/window_len]);           
[Err_index, force_estimate_index] = e_sifir_tst(x_plus_index, feature_test_norm, force_test_index, [Q D Tol ii], [2/window_len 2/window_len]);               
[Err_middle, force_estimate_middle] = e_sifir_tst(x_plus_middle, feature_test_norm, force_test_middle, [Q D Tol ii], [2/window_len 2/window_len]);               
[Err_ring, force_estimate_ring] = e_sifir_tst(x_plus_ring, feature_test_norm, force_test_ring, [Q D Tol ii], [2/window_len 2/window_len]);               
[Err_little, force_estimate_little] = e_sifir_tst(x_plus_little, feature_test_norm, force_test_little, [Q D Tol ii], [2/window_len 2/window_len]);               
            
