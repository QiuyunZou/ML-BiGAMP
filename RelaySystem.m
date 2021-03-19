function obj=RelaySystem(Input)

%% load parameters
Kp=Input.K(1);
Kd=Input.K(2);
N=Input.N;
M=Input.M;
P=Input.P;
nuw1=Input.nuw1;
nuw2=Input.nuw2;
n_bit=Input.n_bit;
ADC_switch=Input.ADC_switch;

%% Generate x
[X,xo,informationBit]=Gen_QAM(Input, N, Kp+Kd);

X_pilot=X(:,1:Kp);
X_Data=X(:,end-Kd+1:end);

%% System
H=(randn(M,N)+1j*randn(M,N))/sqrt(2*N);                 %Channel of first hop
C=(randn(P,M)+1j*randn(P,M))/sqrt(2*M);                 %channel of second hop
W1=sqrt(nuw1/2)*(randn(M,Kp+Kd)+1j*randn(M,Kp+Kd));      %Gaussian noise
W2=sqrt(nuw2/2)*(randn(P,Kp+Kd)+1j*randn(P,Kp+Kd));      %Gaussian noise
S=H*X+W1;
Y=C*S+W2;
[tilde_Y,quan_step]=Quantization(Y,n_bit,ADC_switch);


%% load parameters
obj.X=X;
obj.xo=xo;
obj.tilde_Y=tilde_Y;
obj.Y=Y;
obj.quan_step=quan_step;
obj.H=H;
obj.C=C;
obj.informationBit=informationBit;
obj.X_pilot=X_pilot;
obj.X_Data=X_Data;

end