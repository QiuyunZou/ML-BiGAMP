clc;
clear all;


%% Parameters
K=1000;
N1=10:5:100;
N2=1000;
N3=500;
epsilon=0.01:0.01:0.3;
 

TestNum=1e1;      % Test number
IterNum=150;      % Iteration number;
SNR=50;




%% load parameters
Input.K=K;
Input.N2=N2;
Input.N3=N3;
Input.IterNum=IterNum;


%% Run
parfor_progress(TestNum*length(N1)*length(epsilon));
Testz=zeros(length(N1),length(epsilon));
Input.nuw1=10^(-SNR/10); 
Input.nuw2=10^(-SNR/10);
for a1=1:length(N1)
    Input.N1=N1(a1);
    for a2=1:length(epsilon)
        Input.epsilon=epsilon(a2);
        Detez=zeros(TestNum,1);
        parfor ii=1:TestNum
            obj=BiSystem(Input); 
            Detez(ii,1)=MLBiGAMP_MC(Input,obj);
            parfor_progress;   
        end
        Testz(a1, a2)=mean(Detez);
    end
end
parfor_progress(0);    


Iter=1:IterNum;

figure 
pcolor(epsilon, N1, Testz);
xlabel('Sampling ratio \epsilon');
ylabel('Rank N1');
grid on
colormap gray
colorbar






