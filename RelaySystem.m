function obj=RelaySystem(Input)

%% load parameters
Kp=Input.K(1);
Kd=Input.K(2);
N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
nuw1=Input.nuw1;
nuw2=Input.nuw2;
n_bit=Input.n_bit;
ADC_switch=Input.ADC_switch;

%% Generate x
[X,xo,informationBit]=Gen_QAM(Input, N1, Kp+Kd);
X_pilot=X(:,1:Kp);
X_Data=X(:,end-Kd+1:end);

%% System
H1=(randn(N2,N1)+1j*randn(N2,N1))/sqrt(2*N1);                %Channel of first hop
H2=(randn(N3,N2)+1j*randn(N3,N2))/sqrt(2*N2);                %channel of second hop
W1=sqrt(nuw1/2)*(randn(N2,Kp+Kd)+1j*randn(N2,Kp+Kd));      %Gaussian noise
W2=sqrt(nuw2/2)*(randn(N3,Kp+Kd)+1j*randn(N3,Kp+Kd));      %Gaussian noise
X2=H1*X+W1;
Y=H2*X2+W2;
[tilde_Y,quan_step]=Quantization(Y,n_bit,ADC_switch);

obj.X=X;
obj.xo=xo;
obj.tilde_Y=tilde_Y;
obj.Y=Y;
obj.quan_step=quan_step;
obj.H1=H1;
obj.H2=H2;
obj.informationBit=informationBit;
obj.X_pilot=X_pilot;
obj.X_Data=X_Data;

end