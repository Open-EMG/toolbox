clear all;
close all;

% representative: subject 1, session 2, middle finger, second trial 

subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

path_data=['D:\jiangxinyu\Open Access HDsEMG DataSet Analysis\hyser_dataset']; % change path_data to the location of the dataset
path_results='..\results\decomposition results\1DoF'; % change path_results to the location to save decomposition results

window_len=0.02; % 0.02 s sliding window to extract features
step_len=0.02; % 0.02 s step length of the sliding window to extract features
f_cutoff=10; % filter force data using a low pass filter with 10 Hz cut off frequency
fs_emg=2048;
fs_force=100;
R = 4; % Extension parameter
M = 300; % FastICA iteration number

threshold = 0.8; % threshold to select motor units

for i=1:20
    for session=1:2        
        
        force_finger_1dof=load_1dof(path_data,subject{1,i},session,'force');
        emg_finger_1dof=load_1dof(path_data,subject{1,i},session,'raw');
        
        mvc=get_mvc(path_data,subject{1,i},num2str(session));
        force_norm=normalize_force(force_finger_1dof,mvc);
        force_norm_preprocess=preprocess_force(force_norm,window_len,step_len,f_cutoff,fs_force,fs_emg);
        
        for u=1:5
            for v=1:3        
                
                SpikeTrainGoodSelect=[];
                SILSelect=[];
                for muscle_idx=1:2
                    i
                    session
                    u
                    v
                    muscle_idx
                    emg=emg_finger_1dof{u,v}(:,(muscle_idx-1)*128+1:muscle_idx*128);
                    
                    load([path_results,'/B_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat']);
                    load([path_results,'/SpikeTrain_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat']);
                    load([path_results,'/GoodIndex_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat']);
                    load([path_results,'/SpikeTrainGood_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat']);
                    load([path_results,'/sGood_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat']);
                    load([path_results,'/SIL_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat']);
                    
                    idx=find((SIL>threshold) & (SIL<0.99));
                    SpikeTrainGoodSelect=[SpikeTrainGoodSelect,SpikeTrainGood(:,idx)];
                    SILSelect=[SILSelect,SIL(idx)];
                    if(muscle_idx==1)
                        NumMU_muscle1(u,v,i,session)=size(SpikeTrainGood(:,idx),2);
                    else
                        NumMU_muscle2(u,v,i,session)=size(SpikeTrainGood(:,idx),2);
                    end
                end
                NumMU(u,v,i,session)=size(SpikeTrainGoodSelect,2);
                SILMean(u,v,i,session)=mean(SILSelect);
                if(size(SpikeTrainGoodSelect,2)~=0)
                    [FiringRate,TimeWin] = FiringRateEst(SpikeTrainGoodSelect,[1:NumMU(u,v,i,session)]',window_len,step_len,fs_emg,'On');
                    [b,a]= butter(8,10/(1/(TimeWin(2)-TimeWin(1))/2),'low'); 
                    FiringRate = filtfilt(b,a,double(FiringRate));
                    parameter = regress(force_norm_preprocess{u,v}(51:end-51,u),FiringRate(51:end-51,:));% remove the first and last 1 s signal
                    force_regress=FiringRate*parameter;
                    corr_mat=corrcoef([force_regress(51:end-51),force_norm_preprocess{u,v}(51:end-51,u)]);
                    correlation(u,v,i,session)=corr_mat(1,2);
                else
                    correlation(u,v,i,session)=NaN;
                end
            end
        end
    end
end
SIL_subject_mean=squeeze(mean(mean(mean(SILMean,1,'omitnan'),2,'omitnan'),4,'omitnan'));
correlation_subject_mean=squeeze(mean(mean(mean(correlation,1,'omitnan'),2,'omitnan'),4,'omitnan'));
NumMU_muscle1_subject_mean=squeeze(mean(mean(mean(NumMU_muscle1,1,'omitnan'),2,'omitnan'),4,'omitnan'));
NumMU_muscle2_subject_mean=squeeze(mean(mean(mean(NumMU_muscle2,1,'omitnan'),2,'omitnan'),4,'omitnan'));    
