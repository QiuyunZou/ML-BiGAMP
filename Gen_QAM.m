function [X,xo,informationBit]=Gen_QAM(Input,N,K)

mode_size=Input.mode_size;
sym = modem.qammod(2^mode_size);
normal_scal = 1/sqrt((2/3)*(2^mode_size-1)) ;      %QAM normalization 
sym.input='bit';                                   %the type of input data shoud be 'bit'
sym.symbolorder='gray';
informationBit = round(rand(N*K*mode_size,1)) ;    %产生串行二进制序列
informationSym = modulate(sym, informationBit);    %对二进制序列进行QAM调制
X= normal_scal * reshape(informationSym,N,K) ;     %得到归一化的QAM基带信号
xo=Gen_Constellation(mode_size);

end