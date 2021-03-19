function [tilde_Y,quan_step]=Quantization(Y,n_bit,AGC_switch)
if AGC_switch==0
    tilde_Y=Y;
    quan_step=2;
elseif AGC_switch==1
    [M,K]=size(Y);
    Y=reshape(Y,M*K,1);
    Y_real=real(Y);
    Y_imag=imag(Y);
    
    switch n_bit
        case 1
            quan_step=2;
        case 2
             quan_step=sqrt(0.25);
        case 3
            quan_step=sqrt(0.25);
        otherwise
            quan_step=min(max(abs(Y_real)),max(abs(Y_imag)))/2^(n_bit-1);    
    end
    
    DeltaTh = (0:1:2^(n_bit-1)-1) * quan_step ;
    Q_Out = (1:2:2^n_bit-1)*quan_step/2 ;
    
    Y_R = Q_Out(end) * ( abs(Y_real) > DeltaTh(end) ) ;
    for bIdx = length(DeltaTh):-1:2
        Y_R = Y_R + Q_Out(bIdx-1) * ( (abs(Y_real) < DeltaTh(bIdx)) & (abs(Y_real) >= DeltaTh(bIdx-1)) ) ;
    end
    Y_R = sign(Y_real) .* Y_R ;
    
    Y_I = Q_Out(end) * ( abs(Y_imag) > DeltaTh(end) ) ;
    for bIdx = length(DeltaTh):-1:2
        Y_I = Y_I + Q_Out(bIdx-1) * ( (abs(Y_imag) < DeltaTh(bIdx)) & (abs(Y_imag) >= DeltaTh(bIdx-1)) ) ;
    end
    Y_I = sign(Y_imag) .* Y_I ;
    
    tilde_Y = Y_R + 1j*Y_I ;   
    tilde_Y=reshape(tilde_Y,M,K);
else
    disp('error');
end
end