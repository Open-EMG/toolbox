function force_preprocessed=preprocess_force(force,window_len,step_len,f_cutoff,fs_force,fs_emg)

[row_num,column_num]=size(force);

for i=1:row_num
    for j=1:column_num
        force_tmp=force{i,j};
        [b,a]= butter(8,10/(fs_force/2),'low'); 
        force_tmp_filter = filtfilt(b,a,double(force_tmp));
        force_tmp_filter_resample = resample(force_tmp_filter,fs_emg,fs_force);
        [Nsample,Nchannel]=size(force_tmp_filter_resample);
        window_sample=floor(window_len*fs_emg);
        step_sample=floor(step_len*fs_emg);
        idx=0;
        for u=1:step_sample:(Nsample-window_sample+1)
            idx=idx+1;
            force_preprocessed_tmp(idx,:)=mean(force_tmp_filter_resample(u:u+window_sample-1,:));
        end
        force_preprocessed{i,j}=force_preprocessed_tmp;
    end
end