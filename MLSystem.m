function obj=MLSystem(Input)

%% load parameters
N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
N4=Input.N4;
nuw=Input.nuw;
n_bit=Input.n_bit;
ADC_switch=Input.ADC_switch;

%% Generate x
[X,xo,informationBit]=Gen_QAM(Input, N1, 1);

%% System
H1=(randn(N2, N1)+1j*randn(N2, N1))/sqrt(2*N1);    
H2=(randn(N3, N2)+1j*randn(N3, N2))/sqrt(2*N2); 
H3=(randn(N4, N3)+1j*randn(N4, N3))/sqrt(2*N3); 
              
W1=sqrt(nuw/2)*(randn(N2, 1)+1j*randn(N2, 1));        
W2=sqrt(nuw/2)*(randn(N3, 1)+1j*randn(N3, 1));
W3=sqrt(nuw/2)*(randn(N4, 1)+1j*randn(N4, 1));
       
X2=H1*X+W1;
X3=H2*X2+W2;
Y=H3*X3+W3;

[tilde_Y, quan_step, DeltaTh]=Quantization(Y, n_bit, ADC_switch);

obj.X=X;
obj.xo=xo;
obj.tilde_Y=tilde_Y;
obj.quan_step=quan_step;
obj.H1=H1;
obj.H2=H2;
obj.H3=H3;
obj.informationBit=informationBit;
obj.DeltaTh=DeltaTh;

end