function [BER, SER] = BER_Calculation(X_Est, obj, Input) 


informationBit=obj.informationBit;
modType=Input.modType;
mod_size=Input.mod_size;

[N1,L2] = size(X_Est);

X_Est = reshape(X_Est, N1*L2,1) ; 


switch modType
    case 'QAM',
        normal_scal = 1/sqrt((2/3)*(2^mod_size-1)) ;   % normalization scale: QPSK = 1/sqrt(2) 16QAM = 1/sqrt(10)
        sym = modem.qamdemod(2^mod_size);
    case 'PSK',
        normal_scal = 1 ;
        sym = modem.pskdemod(2^mod_size);
end

% normal_scal = 1/sqrt(mean(abs(qammod(0:2^mod_size-1,2^mod_size)).^2));   % normalization scale: QPSK = 1/sqrt(2) 16QAM = 1/sqrt(10)
X_Est = (1/normal_scal) * X_Est;
 
% sym = modem.qamdemod(2^mod_size);
sym.outputtype = 'bit';
sym.symbolorder = 'gray';
sym.decisiontype = 'hard decision';
z = demodulate(sym, X_Est);
n_bits = N1*L2*mod_size ; 
BER = (n_bits-sum(informationBit==z))/n_bits ;

%% SER
switch modType
    case 'QAM', 
        sym = modem.qammod(2^mod_size);
    case 'PSK',  
        sym = modem.pskmod(2^mod_size);
end 
sym.inputtype = 'bit';
sym.symbolorder = 'gray';
Est_Sym = modulate(sym, z);
informationSym = modulate(sym, informationBit);
n_Syms = N1*L2 ; 
SER = length(find(informationSym-Est_Sym))/n_Syms ;

end
