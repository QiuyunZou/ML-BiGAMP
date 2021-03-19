function [hatx, MSEx, MSEh]=PilotOnly(Input,obj)

%% Channel estimation 
   [hath, MSEh]=FirstPhase(Input, obj);
   
%% Data detection
   [hatx, MSEx]=SecondPhase(Input, obj, hath);

end