function [X,xo,informationBit]=Gen_QAM(Input,N,K)

mode_size=Input.mode_size;
sym = modem.qammod(2^mode_size);
normal_scal = 1/sqrt((2/3)*(2^mode_size-1)) ;      %QAM normalization 
sym.input='bit';                                   %the type of input data shoud be 'bit'
sym.symbolorder='gray';
informationBit = round(rand(N*K*mode_size,1)) ;    %�������ж���������
informationSym = modulate(sym, informationBit);    %�Զ��������н���QAM����
X= normal_scal * reshape(informationSym,N,K) ;     %�õ���һ����QAM�����ź�
xo=Gen_Constellation(mode_size);

end