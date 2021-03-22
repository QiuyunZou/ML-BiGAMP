function [tilde_y, quan_step, DeltaTh]=Quantization(y,n_bit,ADC_switch)
if ADC_switch==0
    tilde_y=y;
    quan_step=0;
    DeltaTh=0;
elseif ADC_switch==1
    y_real=real(y);
    y_imag=imag(y);
    
    if n_bit==1
        quan_step=2;
    elseif n_bit==2 
        quan_step=1/2;
    else
        quan_step=1/(2^(n_bit-1));    
    end
    
    DeltaTh = (0:1:2^(n_bit-1)-1) * quan_step ;    %ADC输入的电平判断
    Q_Out = (1:2:2^n_bit-1)*quan_step/2 ;          %ADC输出电平
    
    y_R = Q_Out(end) * ( abs(y_real) >=DeltaTh(end) ) ;
    for bIdx = length(DeltaTh):-1:2
        y_R = y_R + Q_Out(bIdx-1) * ( (abs(y_real) < DeltaTh(bIdx)) & (abs(y_real) >= DeltaTh(bIdx-1)) ) ;
    end
    y_R = sign(y_real) .* y_R ;
    
    y_I = Q_Out(end) * ( abs(y_imag) >= DeltaTh(end) ) ;
    for bIdx = length(DeltaTh):-1:2
        y_I = y_I + Q_Out(bIdx-1) * ( (abs(y_imag) < DeltaTh(bIdx)) & (abs(y_imag) >= DeltaTh(bIdx-1)) ) ;
    end
    y_I = sign(y_imag) .* y_I ;
    
    tilde_y = y_R + 1j*y_I ;    
else
    disp('error');
end
end