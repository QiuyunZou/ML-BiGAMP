function   obj=MIMO_system(Input)

mod_size=Input.mod_size;
N1=Input.N1;
N2=Input.N2;
bit=Input.bit;
nuw1=Input.nuw1;
ADC_switch=Input.ADC_switch;

%% Generate xo
obj.xo=Gen_Constellation(mod_size); % Generate original constellation sets

%% Generate x
sym = modem.qammod(2^mod_size);
normal_scal = 1/sqrt((2/3)*(2^mod_size-1)) ;      %QAM normalization 
sym.input='bit';                                  %the type of input data shoud be 'bit'
sym.symbolorder='gray';
informationBit = round(rand(N1*mod_size,1)) ;      %产生串行二进制序列
informationSym = modulate(sym, informationBit);   %对二进制序列进行QAM调制
x = normal_scal * reshape(informationSym,N1,1) ;   %得到归一化的QAM基带信号

%% Channel
H=(randn(N2,N1)+1j*randn(N2,N1))/sqrt(2*N1);

%% Noise
w=sqrt(nuw1/2)*(randn(N2,1)+1j*randn(N2,1));   %产生高斯噪声

%% Uncoded system
z=H*x;
y=z+w;
[y,quan_step,DeltaTh]=Quantization(y,bit,ADC_switch);

%% load parameters
obj.x=x;
obj.H=H;
obj.z=z;
obj.y=y;
obj.informationBit=informationBit;
obj.quan_step=quan_step;
obj.DeltaTh=DeltaTh;

end