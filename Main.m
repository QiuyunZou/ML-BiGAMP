clc;
clear all;

%% Parameters setting 
N=50;
M=150;
P=300;
K_tem=500;
beta_p=[0.1:0.004:0.2, 0.25:0.05:0.95,0.99];

n_bit=3;
mes=0.7;          % damping factor;
mode_type='QAM';      
mode_size=2;      % QPSK
TestNum=1e1;      % Test number
IterNum=20;       % Iteration number;
ADC_switch=1;     % 1 on; 0 off.
SNR_dB=5;   


%% load parameters
Input.N=N;
Input.M=M;
Input.P=P; 
Input.mes=mes;
Input.n_bit=n_bit;
Input.IterNum=IterNum;
Input.mode_type=mode_type;
Input.mode_size=mode_size;
Input.ADC_switch=ADC_switch;
Input.nuw1=10^(-SNR_dB/10);
Input.nuw2=10^(-SNR_dB/10);


%% Run
tic;
parfor_progress(TestNum*length(beta_p));
for Index=1:length(beta_p)
    Kp=round(beta_p(Index)*K_tem);
    Kd=K_tem-Kp;
    Input.K=[Kp, Kd];
    parfor ii=1:TestNum
        obj=RelaySystem(Input); 
        [Pilot_x, Pilot_MSE_x, Pilot_h]=PilotOnly(Input,obj);
        [JCD_x, JCD_MSE_x, JCD_h ]=JCD(Input,obj);
        
        [Pilot_BER(Index,ii), Pilot_SER(Index ,ii)] = BER_Calculation(Pilot_x, Input, obj);
        [JCD_BER(Index,ii), JCD_SER(Index ,ii)] = BER_Calculation(JCD_x, Input, obj);
        
        Pilot_H(Index, ii)=Pilot_h(end);
        JCD_H(Index, ii)=JCD_h(end);
        parfor_progress;         % Count
    end
    Pilot_BER_Mean(Index,1)=mean(Pilot_BER(Index,:));
    JCD_BER_Mean(Index,1)=mean(JCD_BER(Index,:));
      
    Pilot_Mean_H(Index, 1)=mean(Pilot_H(Index,:));
    JCD_Mean_H(Index, 1)=mean(JCD_H(Index,:));
end
parfor_progress(0);       % Clean count
toc;




figure(1)
semilogy(beta_p,  Pilot_Mean_H, '-.r');  hold on;
semilogy(beta_p,  JCD_Mean_H, '--b');  hold on;
legend('Pilot-only', 'JCD'); hold on;
xlabel('beta_p=Kp/K'); hold on;
ylabel('MSE'); hold on ;






