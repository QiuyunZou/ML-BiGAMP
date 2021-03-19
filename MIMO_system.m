function   obj=MIMO_system(Input)

mod_size=Input.mod_size;
N=Input.N;
M=Input.M;
bit=Input.bit;
nuw=Input.nuw;
AGC_switch=Input.AGC_switch;

%% Generate xo
obj.xo=Gen_Constellation(mod_size); % Generate original constellation sets

%% Generate x
sym = modem.qammod(2^mod_size);
normal_scal = 1/sqrt((2/3)*(2^mod_size-1)) ;      %QAM normalization 
sym.input='bit';                                  %the type of input data shoud be 'bit'
sym.symbolorder='gray';
informationBit = round(rand(N*mod_size,1)) ;      %产生串行二进制序列
informationSym = modulate(sym, informationBit);   %对二进制序列进行QAM调制
x = normal_scal * reshape(informationSym,N,1) ;   %得到归一化的QAM基带信号

%% Channel
H=(randn(M,N)+1j*randn(M,N))/sqrt(2*N);

%% Noise
w=sqrt(nuw/2)*(randn(M,1)+1j*randn(M,1));   %产生高斯噪声

%% Uncoded system
y=H*x+w;
[tilde_y,quan_step]=Quantization(y,bit,AGC_switch);

%% load parameters
obj.x=x;
obj.H=H;
obj.y=tilde_y;
obj.informationBit=informationBit;
obj.quan_step=quan_step;

end