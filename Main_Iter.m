clc;
clear all;

%% Parameters Setting
N=512;
M=512;
mes=0.95;
AGC_switch=1;
TestNum=1e2;
IterNum=20;
modType='QAM';
mod_size=2;
bit=1;
snr=12;

%% Load Parameters
Input.N=N;
Input.M=M;
Input.AGC_switch=AGC_switch;
Input.IterNum=IterNum;
Input.modType=modType;
Input.mod_size=mod_size;
Input.bit=bit;
Input.mes=mes;
Input.nuw=10^(-snr/10);


parfor_progress(TestNum);
for Index=1:TestNum
    obj=MIMO_system(Input);
    GAMP_MSE(:,Index)=MLBiGAMP(obj,Input);
    parfor_progress;   
end
parfor_progress(0); 

for kk=1:IterNum
    GAMP_MSE_mean(kk)=mean(GAMP_MSE(kk,:));
end


MSE =MLBiGAMP_SE(Input,obj);

Dis=1:IterNum;
semilogy(Dis,  GAMP_MSE_mean, '-*b');   hold on;
semilogy(Dis,  MSE, 'or');   hold on;
legend('GAMP'); 
hold on

xlabel('SNR');
ylabel('BER');
