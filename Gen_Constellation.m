function xo=Gen_Constellation(mod_size)

normal_scal = 1/sqrt((2/3)*(2^mod_size-1)) ;      %Normalization factor
sym = modem.qammod(2^mod_size);
sym.inputtype = 'integer';   %��������Ϊ��������
sym.symbolorder = 'gray';    %������
xo = normal_scal * modulate(sym, 0:2^mod_size-1); %�õ������㼯��

end