clear all;
close all;

path_data='D:\jiangxinyu\Open Access HDsEMG DataSet Analysis\hyser_dataset'; % change path_data to the location of the dataset
path_results='..\results\decomposition results\1DoF'; % change path_results to the location to save decomposition results

fs_emg=2048;
R = 4; % Extension parameter
M = 300; % FastICA iteration number

threshold = 0.6; % threshold to select motor units

subject='01';
session=1;
finger=2;%index finger
sample=1; % take the sample acquired in the first trial as an example;

emg_finger_1dof_all=load_1dof(path_data,subject,session,'raw');
emg_finger_1dof=emg_finger_1dof_all{finger,sample};

muscle_idx=1; % take extensor muscle as an example

emg=emg_finger_1dof(:,(muscle_idx-1)*128+1:muscle_idx*128);
[emg_extend,W] = SimEMGProcessing(emg,'R',R,'WhitenFlag','On','SNR','Inf');  % WhitenFlag: whiten the signal or not; SNR: add XdB noises or not; Inf: do not add noises
[s,B,SpikeTraintemp] = FastICA(emg_extend,M,fs_emg);
[SpikeTrain,GoodIndex] = MUReplicasRemoval(SpikeTraintemp,s,fs_emg); 
SpikeTrainGood = SpikeTrain(:,GoodIndex);
sGood=s(:,GoodIndex);
SIL = SILCal(sGood,fs_emg); %calculate the SIL value of each source 
NumMU = length(SIL((SIL>threshold) & (SIL<0.99)));
MeanSIL = mean(SIL((SIL>threshold) & (SIL<0.99)));
