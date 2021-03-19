clc;
clear all;

% This is a program of signal reconstruction of two-hop wireless communication: X2=H1*[Xp, Xd]+W1, Y=Qc(H2*X2+W2)
% where Xp is the pilot and Xd is the user data. Both Xp and Xd is QPSK symbols. 
% Xp:Np*K, Xd: Nd*K; H1: N2*N1; X2: N2*K; H2: N3*N2; Y: N3*K.
% H2 is knwon, H1 is to be estimated. 


%% Parameters
Kp=100;           % length of pilot
Kd=500-Kp;        % length of data
K=[Kp,Kd];    
N1=50;
N2=200;
N3=400;

n_bit=3;          % the number of ADC bit
mes=0.75;         % damping factor;    
mode_size=2;      % QPSK
TestNum=1e4;      % Test number
IterNum=30;       % Iteration number;
ADC_switch=0;     % 0: infinite bit; 1: finite bits.
SNR1=8;           % SNR of 1-st layer;  
SNR2=8;           % SNR of 2-ed layer;


%% load parameters
Input.K=K;
Input.N1=N1;
Input.N2=N2;
Input.N3=N3; 
Input.mes=mes;
Input.n_bit=n_bit;
Input.IterNum=IterNum;
Input.mode_size=mode_size;
Input.ADC_switch=ADC_switch;
Input.nuw1=10^(-SNR1/10); 
Input.nuw2=10^(-SNR2/10); 


%% Run
tic;
parfor_progress(TestNum);   % Counter
parfor ii=1:TestNum
    obj=RelaySystem(Input); 
    [MSE_x(:,ii), MSE_h(:,ii),q_h(ii)]=MLBiGAMP(Input,obj);
    parfor_progress;   
end
parfor_progress(0);         % clear 

for ii=1:IterNum
    MSE_x_mean(ii,1)=mean(MSE_x(ii,:));
    MSE_h_mean(ii,1)=mean(MSE_h(ii,:));
end

[MSE_x,MSE_h]=MLBiGAMP_SE(Input,q_h);
toc;

Iter=1:IterNum;

figure(1)
plot(Iter,10*log10(MSE_x_mean), '-ob');  hold on;
plot(Iter,10*log10(MSE_x), '*r');        hold on;
legend('ML-BiGAMP','SE');  hold on;
xlabel('Iteration');       hold on;
ylabel('MSE');             hold on;





