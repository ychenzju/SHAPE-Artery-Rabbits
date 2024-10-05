function [set_num] = count_sets(arrX)

set_num = 0;
set_start = 0;
set_stop = 0;
for i = 1:length(arrX)
    if arrX(i)>0 && set_start==0 && set_stop==0
        set_start = 1;
        set_num = set_num + 1;
    elseif arrX(i)>0 && set_start==1 && set_stop==0
        continue
    elseif arrX(i)>0 && set_start==0 && set_stop==1
        set_stop = 0;
    elseif arrX(i)>0 && set_start==1 && set_stop==1
        set_start = 0;
        set_stop = 0;
    elseif arrX(i)<=0 && set_start==0 && set_stop==0
        continue
    elseif arrX(i)<=0 && set_start==1 && set_stop==0
        set_stop = 1;
    elseif arrX(i)<=0 && set_start==0 && set_stop==1
        set_stop = 0;
    elseif arrX(i)<=0 && set_start==1 && set_stop==1
        set_start = 0;
        set_stop = 0;
    end
    
    if set_start==1 && set_stop==1
        set_start= 0;
        set_stop = 0;
    end
end
        