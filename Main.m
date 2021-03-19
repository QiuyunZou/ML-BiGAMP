clc;
clear all;

%% Parameters setting 
Kp=100;
Kd=400;
K=[Kp Kd];
N=50;
M=150;
P=300;

SNR_dB=0:1:20;

n_bit=1;
mes=0.7;          % damping factor;
mode_type='QAM';      
mode_size=2;      % QPSK
TestNum=1e1;      % Test number
IterNum=30;       % Iteration number;
ADC_switch=1;     % 1 on; 0 off. 


%% load parameters
Input.K=K;
Input.N=N;
Input.M=M;
Input.P=P; 
Input.mes=mes;
Input.n_bit=n_bit;
Input.IterNum=IterNum;
Input.mode_type=mode_type;
Input.mode_size=mode_size;
Input.ADC_switch=ADC_switch;



%% Run
tic;
parfor_progress(TestNum*length(SNR_dB));
for Index=1:length(SNR_dB)
    Input.nuw1=10^(-SNR_dB(Index)/10);
    Input.nuw2=10^(-SNR_dB(Index)/10);
    parfor ii=1:TestNum
        obj=RelaySystem(Input); 
        [Pilot_x, Pilot_MSEx, Pilot_h]=PilotOnly(Input,obj);
        [JCD_x, JCD_MSEx, JCD_h ]=JCD(Input,obj);
        [PerfectCSI_x, PerfectCSI_MSEx]=PerfectCSI(Input,obj);
                
        [Pilot_BER(Index,ii), Pilot_SER(Index ,ii)] = BER_Calculation(Pilot_x, Input, obj);
        [JCD_BER(Index,ii), JCD_SER(Index ,ii)] = BER_Calculation(JCD_x, Input, obj);
        [PerfectCSI_BER(Index,ii), PerfectCSI_SER(Index ,ii)] = BER_Calculation(PerfectCSI_x, Input, obj);
        parfor_progress;         % Count
    end
    Pilot_BER_Mean(Index,1)=mean(Pilot_BER(Index,:));
    JCD_BER_Mean(Index,1)=mean(JCD_BER(Index,:));
    PerfectCSI_BER_Mean(Index,1)=mean(PerfectCSI_BER(Index,:));
end
parfor_progress(0);       % Clean count
toc;



figure(1)
semilogy(SNR_dB,  Pilot_BER_Mean, '-*b');  hold on;
semilogy(SNR_dB,  JCD_BER_Mean, '-or');  hold on;
semilogy(SNR_dB,  PerfectCSI_BER_Mean, '-sk');  hold on;
legend('Pilot-only', 'JCD', 'Perfect-CSI'); hold on;
xlabel('SNR[dB]'); hold on;
ylabel('BER'); hold on ;







