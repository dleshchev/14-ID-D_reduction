function [ time_str ] = time_num2str( time_num )
%TIME_NUM2STR Summary of this function goes here
%   Detailed explanation goes here

A = log10(abs(time_num));

if      A>=-12 && A<-9
    time_str = [num2str(time_num*1e12) ' ps'];
elseif  A>=-9 && A<-6
     time_str = [num2str(time_num*1e9) ' ns'];
elseif  A>=-6 && A<-3
     time_str = [num2str(time_num*1e6) ' us'];
elseif  A>=-3 && A<0
     time_str = [num2str(time_num*1e3) ' ms'];
elseif A>=0
     time_str = [num2str(time_num) ' s'];
end


end

