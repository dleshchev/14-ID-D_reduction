function [ time_num ] = time_str2num( time_str )
%TIME_STR2NUM Summary of this function goes here
%   Detailed explanation goes here

if     strfind(time_str,'fs')
    time_num = get_time_num(time_str,'fs')*1e-15;

elseif strfind(time_str,'ps')
    time_num = get_time_num(time_str,'ps')*1e-12;

elseif strfind(time_str,'ns')
    time_num = get_time_num(time_str,'ns')*1e-9;

elseif strfind(time_str,'us')
    time_num = get_time_num(time_str,'us')*1e-6;

elseif strfind(time_str,'ms')
    time_num = get_time_num(time_str,'ms')*1e-3;

elseif strfind(time_str,'s')
    time_num = get_time_num(time_str,'s');
else
    time_num = str2double(time_str);
end   
    


    function time_num = get_time_num(time_str,unit_str)
        idx = strfind(time_str,unit_str);
        time_num = str2double(time_str(1:idx-1));
    end
end

