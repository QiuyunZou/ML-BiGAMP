clc;
clear all;

%% Parameters Setting
N=1024;
M=512;
mes=0.95;
ADC_switch=1;
TestNum=1e2;
IterNum=20;
bit=1;
snr=10;
rho=0.05;

%% Load Parameters
Input.N=N;
Input.M=M;
Input.ADC_switch=ADC_switch;
Input.IterNum=IterNum;
Input.rho=rho;
Input.bit=bit;
Input.mes=mes;
Input.nuw=10^(-snr/10);

for Index=1:TestNum
    obj=MIMO_system(Input);
    MLBiGAMP_MSE(:,Index)=MLBiGAMP(obj,Input);
    if (mod(Index,TestNum/10)==0)
        disp(Index/TestNum*10);
    end
end

for kk=1:IterNum
    GAMP_MSE_mean(kk)=mean(MLBiGAMP_MSE(kk,:));
end


MSE =GAMP_SE(Input,obj);

Dis=1:IterNum;
plot(Dis,  10*log10(GAMP_MSE_mean), '-*b');   hold on;
plot(Dis,  10*log10(MSE), 'or');   hold on;
legend('GAMP', 'SE'); 
hold on

xlabel('SNR');
ylabel('MSE(dB)');
