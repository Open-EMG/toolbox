clear all;
close all;

subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};

path_data='D:\jiangxinyu\Open Access HDsEMG DataSet Analysis\hyser_dataset'; % change path_data to the location of the dataset
path_results='..\results\decomposition results\1DoF'; % change path_results to the location to save decomposition results
mkdir(path_results);
          
fs_emg=2048;
R = 4; % Extension parameter
M = 300; % FastICA iteration number

threshold = 0.6; % shreshold to select motor units

parfor i=1:20
    for session=1:2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        emg_finger_1dof=load_1dof(path_data,subject{1,i},session,'raw');
        for u=1:5 %%%%%%%%%%%%%%%%%% u=1:5
            for v=1:3
                for muscle_idx=1:2
                    i
                    session
                    u
                    v
                    muscle_idx
                    emg=emg_finger_1dof{u,v}(:,(muscle_idx-1)*128+1:muscle_idx*128);
                    [emg_extend,W] = SimEMGProcessing(emg,'R',R,'WhitenFlag','On','SNR','Inf');  % WhitenFlag: whiten the signal or not; SNR: add XdB noises or not; Inf: do not add noises
                    [s,B,SpikeTraintemp] = FastICA(emg_extend,M,fs_emg);
                    [SpikeTrain,GoodIndex] = MUReplicasRemoval(SpikeTraintemp,s,fs_emg); 
                    SpikeTrainGood = SpikeTrain(:,GoodIndex);
                    sGood=s(:,GoodIndex);
                    SIL = SILCal(sGood,fs_emg); %calculate the SIL value of each source 
                    
                    mySaveB([path_results,'/B_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat'],B);
                    mySaveSpikeTrain([path_results,'/SpikeTrain_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat'],SpikeTrain);
                    mySaveGoodIndex([path_results,'/GoodIndex_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat'],GoodIndex);
                    mySaveSpikeTrainGood([path_results,'/SpikeTrainGood_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat'],SpikeTrainGood);
                    mySavesGood([path_results,'/sGood_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat'],sGood);
                    mySaveSIL([path_results,'/SIL_subject',subject{1,i},'_session',num2str(session),'_task',num2str(u),'_trial',num2str(v),'_R',num2str(R),'_M',num2str(M),'_muscle',num2str(muscle_idx),'.mat'],SIL);
                     
                end
            end
        end
    end
end