clc;
clear all;

%% Parameters Setting
N1=200;
N2=400;
N3=600;
N4=800;
mes=0.8;
ADC_switch=1;
TestNum=1e0;    
IterNum=20;         
mode_size=2;
n_bit=1;
snr=15;

%% Load Parameters
Input.N1=N1;
Input.N2=N2;
Input.N3=N3;
Input.N4=N4;
Input.ADC_switch=ADC_switch;
Input.IterNum=IterNum;
Input.mode_size=mode_size;
Input.n_bit=n_bit;
Input.mes=mes;
Input.nuw=10^(-snr/10);

tic;
parfor_progress(TestNum);
parfor kk=1:TestNum
    obj=MLSystem(Input);
    MLBiGAMP_MSE_Error(kk,:)=MLBiGAMP_CSI(Input, obj);
    MLVAMP_MSE_Error(kk,:)=MLVAMP(Input, obj);
    parfor_progress;      % Count
end
parfor_progress(0);       % Clean count
toc;

for ii=1:IterNum
    MLVAMP_MSE_Mean(ii)=mean(MLVAMP_MSE_Error(:,ii));
    MLBiGAMP_MSE_Mean(ii)=mean(MLBiGAMP_MSE_Error(:,ii));
end



iter=1:IterNum;
semilogy(iter,   MLVAMP_MSE_Mean,  '-*r' );   hold on;
semilogy(iter,   MLBiGAMP_MSE_Mean,  '-ok' );   hold on;
legend('ML-VAMP', 'ML-BiGAMP'); hold on;

xlabel('Iter');
ylabel('NMSE');
