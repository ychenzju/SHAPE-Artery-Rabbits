function [tt,ecg,bp,track_time] = getBPfromDB(DB,time_acq,time_start,time_length,time_delta)
% 
% time_begin = '11:58:49'; %'12:10:10'; %'10:50:54';
% time_length = 40;
% DB = data;

eps = 1e-3;
timestrs1 = split(time_acq,':');
timestrs2 = split(time_delta,':');

hour_begin = str2num(timestrs1{1});
minute_begin = str2num(timestrs1{2});
second_begin = str2num(timestrs1{3});

hour_delta = str2num(timestrs2{1});
minute_delta = str2num(timestrs2{2});
second_delta = str2num(timestrs2{3});

if second_begin-second_delta<0
    second_begin = second_begin - second_delta + 60;
    minute_begin = minute_begin -1;
else
    second_begin = second_begin - second_delta;
end

if minute_begin-minute_delta<0
    minute_begin = minute_begin - minute_delta + 60;
    hour_begin = hour_begin -1;
else
    minute_begin = minute_begin - minute_delta;
end
   
if second_begin+time_start>=60
    second_begin = second_begin + time_start - 60;
    minute_begin = minute_begin + 1;
    if minute_begin>=60
        minute_begin = minute_begin - 60;
        hour_begin = hour_begin + 1;
    end
else
    second_begin = second_begin + time_start;
end

list = abs(DB.date_time(:,4)-hour_begin)<eps & ...
    abs(DB.date_time(:,5)-minute_begin)<eps & ...
    abs(DB.date_time(:,6)-second_begin)<eps;

indices = find(list);
index_begin = indices(1);
interval = DB.Interval;

L = floor(time_length/interval);
tt = (0:L-1)' * interval + time_start;
ecg = DB.ECG(index_begin:index_begin+L-1);
bp  = DB.BP(index_begin:index_begin+L-1);
track_time = [num2str(hour_begin),':',num2str(minute_begin),':',num2str(second_begin)];




