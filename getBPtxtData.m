function data = getBPxlsData(filename)
% path_str = 'C:\MyData\XieHe\PT-Trial\All-Rabbits-Data-WY\XR003\';
% filename = [path_str,'XRtest.txt'];

failed = false;
if isfile(filename)
    fileId = fopen(filename,'r');
    ended = 0; row = 0;
    L = 9900000;
    date_time = nan(L,6);
    sectime = nan(L,1);
    ms  = nan(L,1);
    p1 = nan(L,1);
    p2 = nan(L,1);
    data_index = 0;
    next_is_data = 0;
    
    time_begin_ok = 0;
    while ~ended
        
        if mod(data_index, 10000)==1
            disp(data_index)
        end
        
        row = row + 1;
        text = fgetl(fileId);
        if text==-1
            ended = 1;
            continue;
        end
        
        if startsWith(text,'ExcelDateTime') && time_begin_ok==0
            subStrCells = split(text);
            day_time_begin = datevec(datetime([subStrCells{3},' ',subStrCells{4}],...
                'InputFormat','yyyy/MM/dd HH:mm:ss.SSSS'));
            time_begin_ok = 1;
            sectime_begin = 0;
        end
        
        if startsWith(text,'Interval')
            next_is_data = 0;
            subStrCells = split(text);
            Interval = str2double(subStrCells{2});
            continue
        end
        
        if startsWith(text,'BottomValue')
            next_is_data = 1;
            continue
        end
        
        if next_is_data==1
            data_index = data_index + 1;
            subStrCells = split(text);
            second = str2double(subStrCells{1});
            dt = datevec(datetime(subStrCells{2},'Format','yyyy/MM/dd') + seconds(second));
            date_time(data_index,:) = dt;
            ms(data_index) = (second - floor(second)) * 1000;
            p1(data_index) = str2double(subStrCells{3});
            p2(data_index) = str2double(subStrCells{4});
            sectime(data_index) = etime(dt, day_time_begin);
        end
    end
    fclose(fileId);
else
    disp('Incorrect filename');
    failed = true;
end

if failed
    data = struct;
else
    strs = split(filename,'\');
    fn = split(strs{end},'.');
    data.name = fn{1};
    [data.date_time,TF] = rmmissing(date_time);
    data.ECG = p1(~TF); %ECG
    data.BP = p2(~TF); %BP
    data.ms = ms(~TF);
    data.Interval = Interval;
    data.sectime = sectime(~TF);
end