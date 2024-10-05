function datainfo = getCSVformat(filename)
% path_str = 'C:\MyData\XieHe\PT-Trial\2023-05-20-Rabbits\Data\R15\';
% filename = [path_str,'018.csv'];

fileId = fopen(filename,'r');
ended = 0; row = 0;
while ~ended
    text = fgetl(fileId);
    row = row + 1;
    k = strfind(text,'Time');
    if ~isempty(k) && k(1)==1
%         disp('Time is found');
        datainfo.row = row;
        datainfo.time_index = 1;
        datainfo.trace_indices = nan(1,8);
        datainfo.trace_colors = cell(0);
        datainfo.trace_labels = cell(0);
        datainfo.t1_index = nan;
        datainfo.t2_index = nan;
        datainfo.ecg_time_index = nan;
        datainfo.ecg_trace_index = nan;
        C = strsplit(text,',');
        for j = 2:length(C)
            if contains(C{j},'T1')
                datainfo.t1_index = j;
            elseif contains(C{j},'T2')
                datainfo.t2_index = j;
            elseif contains(C{j},'Trace 1')
                datainfo.trace_indices(1) = j; datainfo.trace_colors{1}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{1} = 'Trace#1';
            elseif contains(C{j},'Trace 2')
                datainfo.trace_indices(2) = j; datainfo.trace_colors{2}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{2} = 'Trace#2';
            elseif contains(C{j},'Trace 3')
                datainfo.trace_indices(3) = j; datainfo.trace_colors{3}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{3} = 'Trace#3';
            elseif contains(C{j},'Trace 4')
                datainfo.trace_indices(4) = j; datainfo.trace_colors{4}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{4} = 'Trace#4';
            elseif contains(C{j},'Trace 5')
                datainfo.trace_indices(5) = j; datainfo.trace_colors{5}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{5} = 'Trace#5';
            elseif contains(C{j},'Trace 6')
                datainfo.trace_indices(6) = j; datainfo.trace_colors{6}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{6} = 'Trace#6';
            elseif contains(C{j},'Trace 7')
                datainfo.trace_indices(7) = j; datainfo.trace_colors{7}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{7} = 'Trace#7';
            elseif contains(C{j},'Trace 8')
                datainfo.trace_indices(8) = j; datainfo.trace_colors{8}=extractBetween(C{j},'(',')');
                datainfo.trace_labels{8} = 'Trace#8';
            elseif contains(C{j},'Ecg')
                datainfo.ecg_time_index = j+1;
                datainfo.ecg_trace_index = j+2;
                break;
            end
        end
        ended = 1;
    end
end

if exist('datainfo','var')&& ~isempty(datainfo.trace_indices)
    datainfo.trace_indices = rmmissing(datainfo.trace_indices);
end
fclose(fileId);

