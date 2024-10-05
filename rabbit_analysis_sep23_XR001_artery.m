pathold = path; path(pathold,'.\MatlabFuns');
close all; clc;
clearvars -except data;
format compact
format short g

global plot_on
plot_on =0;


%% define SHAPE export file
base_path = '.\All-Rabbits-Data-WY\';
human_name = 'XR001';
data_path = [base_path,human_name,'\'];
mat_file = [base_path,human_name,'_WY.mat'];
bp_file = [base_path,human_name,'_BP.mat'];

RELOAD_DATA = 'yes'; %'yes' 强制重新读取原始数据。否则从存储的Matlab文件中读取

%% ============== Paper cases ==============
%% Upstream Long-Axis
% [1],[2,3,4,5]
%% Upstream Short-Axis
% [3],[1]
% [4],[1]
% [5],[1]
% [7],[1]
%% Downstream Long-Axis
% [8],[1,2,3,4]
%% Downsteram Short-Axis
% [10],[1]
% [11],[1]
% [12],[1]
% [13],[1]
%% ============== Paper cases ==============

STUDY_CASE_IDS = [8];
STUDY_TRAC_IDS = [1,2,3,4,5];
%% Load all cases
if ~exist(mat_file,'file') || strcmp(RELOAD_DATA,'yes')
    CASES = {};
    default_trace_labels = {};
    
    for caseid = STUDY_CASE_IDS
        start_time_of_interest_default = 0;
        timelength_of_interest_default = 30;
        heart_beat_freq_default = 150/60;
        
        %% define data files
        switch caseid
            %% L2-9 UP
            case 1
                filename = 'L2-9-LX-UP.csv'; hrmin = '10:50:54'; note = 'Long-Axis-UpStream';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'Unknown-1','Proximal-2(Most)','Proximal-3','Proximal-4',...
                    'Distal-5','Distal-6','Distal-7','Distal-8(Most)'};
                visualize_trace_list = [2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 2
                filename = 'L2-9-LX-UP-HRES.csv'; hrmin = '10:52:32'; note = 'Long-Axis-UpStream';
                probe = 'L2-9'; imgmode = 'HRES'; cc = 0.024; mi = 0.16;
                trace_labels = {'Proximal-1(most)','Proximal-2','Proximal-3','Proximal-4',...
                    'Distal-5','Distal-6','Distal-7','Distal-8(Most)'};
                visualize_trace_list = [1,2,3,4,5,6,7,8];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 3
                filename = 'L2-9-SX-UP-01.csv'; hrmin = '10:56:59'; note = 'Short-Axis-UpStream#1';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-UP#1-1','SX-UP#1-2','SX-UP#1-3','SX-UP#1-4',...
                    'SX-UP#1-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 4
                filename = 'L2-9-SX-UP-02.csv'; hrmin = '10:58:14'; note = 'Short-Axis-UpStream#2';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-UP#2-1','SX-UP#2-2','SX-UP#2-3','SX-UP#2-4',...
                    'SX-UP#2-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 5
                filename = 'L2-9-SX-UP-03.csv'; hrmin = '10:59:13'; note = 'Short-Axis-UpStream#3';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-UP#3-1','SX-UP#3-2','SX-UP#3-3','SX-UP#3-4',...
                    'SX-UP#3-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 6
                filename = 'L2-9-SX-UP-04.csv'; hrmin = '11:00:20'; note = 'Short-Axis-UpStream#4';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-UP#4-1','SX-UP#4-2','SX-UP#4-3','SX-UP#4-4',...
                    'SX-UP#4-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 7
                filename = 'L2-9-SX-UP-05.csv'; hrmin = '11:01:23'; note = 'Short-Axis-UpStream#5';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-UP#5-1','SX-UP#5-2','SX-UP#5-3','SX-UP#5-4',...
                    'SX-UP#5-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
                %% L2-9 Down
            case 8
                filename = 'L2-9-LX-DOWN.csv'; hrmin = '11:08:20'; note = 'Long-Axis-DownStream';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'Proximal-1(Most)','Proximal-2','Distal-3','Distal-4(Most)',...
                    'Unknown-5','Unknown-6','Unknown-7'};
                visualize_trace_list = [1,2,3,4];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 9
                filename = 'L2-9-LX-DOWN-HRES.csv'; hrmin = '11:09:16'; note = 'Long-Axis-DownStream';
                probe = 'L2-9'; imgmode = 'HRES'; cc = 0.024; mi = 0.16;
                trace_labels = {'Proximal-1(Most)','Proximal-2','Distal-3','Distal-4(Most)',...
                    'Unknown-5','Unknown-6','Unknown-7'};
                visualize_trace_list = [1,2,3,4];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 10
                filename = 'L2-9-SX-DOWN-01.csv'; hrmin = '11:11:45'; note = 'Short-Axis-DownStream#1';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-DN#1-1','SX-DN#1-2','SX-DN#1-3','SX-DN#1-4',...
                    'SX-DN#1-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 11
                filename = 'L2-9-SX-DOWN-02.csv'; hrmin = '11:12:54'; note = 'Short-Axis-DownStream#2';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-DN#2-1','SX-DN#2-2','SX-DN#2-3','SX-DN#2-4',...
                    'SX-DN#2-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 12
                filename = 'L2-9-SX-DOWN-03.csv'; hrmin = '11:13:50'; note = 'Short-Axis-DownStream#3';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-DN#3-1','SX-DN#3-2','SX-DN#3-3','SX-DN#3-4',...
                    'SX-DN#3-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 13
                filename = 'L2-9-SX-DOWN-04.csv'; hrmin = '11:14:40'; note = 'Short-Axis-DownStream#4';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-DN#4-1','SX-DN#4-2','SX-DN#4-3','SX-DN#4-4',...
                    'SX-DN#4-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 14
                filename = 'L2-9-SX-DOWN-05.csv'; hrmin = '11:15:30'; note = 'Short-Axis-DownStream#5';
                probe = 'L2-9'; imgmode = 'RES1'; cc = 0.024; mi = 0.16;
                trace_labels = {'SX-DN#5-1','SX-DN#5-2','SX-DN#5-3','SX-DN#5-4',...
                    'SX-DN#5-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
                %% C1-6 UP
            case 15
                filename = 'C1-6-LX-UP.csv'; hrmin = '11:58:49'; note = 'Long-Axis-UpStream';
                probe = 'C1-6'; imgmode = 'SHAM'; cc = 0.024; mi = 0.20;
                trace_labels = {'Proximal-1(Most)','Proximal-2','Proximal-3','Proximal-4',...
                    'Distal-5','Distal-6','Distal-7','Distal-8(Most)'};
                visualize_trace_list = [1,2,3,4,5,6,7,8];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 16
                filename = 'C1-6-LX-UP-HRES.csv'; hrmin = '11:59:53'; note = 'Long-Axis-UpStream';
                probe = 'C1-6'; imgmode = 'HRES'; cc = 0.024; mi = 0.20;
                trace_labels = {'Proximal-1(most)','Proximal-2','Proximal-3','Proximal-4',...
                    'Distal-5','Distal-6','Distal-7','Distal-8(Most)'};
                visualize_trace_list = [1,2,3,4,5,6,7,8];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 17
                filename = 'C1-6-SX-UP-01.csv'; hrmin = '12:02:16'; note = 'Short-Axis-UpStream#1';
                probe = 'C1-6'; imgmode = 'SHAM'; cc = 0.024; mi = 0.20;
                trace_labels = {'SX-UP#1-1','SX-UP#1-2','SX-UP#1-3','SX-UP#1-4',...
                    'SX-UP#1-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 18
                filename = 'C1-6-SX-UP-02.csv'; hrmin = '12:03:19'; note = 'Short-Axis-UpStream#2';
                probe = 'C1-6'; imgmode = 'SHAM'; cc = 0.024; mi = 0.20;
                trace_labels = {'SX-UP#2-1','SX-UP#2-2','SX-UP#2-3','SX-UP#2-4',...
                    'SX-UP#2-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 19
                filename = 'C1-6-SX-UP-03.csv'; hrmin = '12:04:14'; note = 'Short-Axis-UpStream#3';
                probe = 'C1-6'; imgmode = 'SHAM'; cc = 0.024; mi = 0.20;
                trace_labels = {'SX-UP#3-1','SX-UP#3-2','SX-UP#3-3','SX-UP#3-4',...
                    'SX-UP#3-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 20
                filename = 'C1-6-SX-UP-04.csv'; hrmin = '12:05:43'; note = 'Short-Axis-UpStream#4';
                probe = 'C1-6'; imgmode = 'SHAM'; cc = 0.024; mi = 0.20;
                trace_labels = {'SX-UP#4-1','SX-UP#4-2','SX-UP#4-3','SX-UP#4-4',...
                    'SX-UP#4-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 21
                filename = 'C1-6-SX-UP-05.csv'; hrmin = '12:06:36'; note = 'Short-Axis-UpStream#5';
                probe = 'C1-6'; imgmode = 'SHAM'; cc = 0.024; mi = 0.20;
                trace_labels = {'SX-UP#5-1','SX-UP#5-2','SX-UP#5-3','SX-UP#5-4',...
                    'SX-UP#5-5','Unknown-6','Unknown-7','Unknown-8'};
                visualize_trace_list = [1,2,3,4,5];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
                %% C1-6 DOWN
            case 22
                filename = 'C1-6-LX-DOWN.csv'; hrmin = '12:10:10'; note = 'Long-Axis-DownStream';
                probe = 'C1-6'; imgmode = 'RES1'; cc = 0.024; mi = 0.20;
                trace_labels = {'Proximal-1(Most)','Proximal-2','Proximal-3','Proximal-4',...
                    'Distal-5','Distal-6','Distal-7','Distal-8(Most)'};
                visualize_trace_list = [1,2,3,4,5,6,7,8];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
            case 23
                filename = 'C1-6-LX-DOWN-HRES.csv'; hrmin = '12:14:14'; note = 'Long-Axis-DownStream';
                probe = 'C1-6'; imgmode = 'HRES'; cc = 0.024; mi = 0.20;
                trace_labels = {'Proximal-1(Most)','Proximal-2','Proximal-3','Proximal-4',...
                    'Distal-5','Distal-6','Distal-7','Distal-8(Most)'};
                visualize_trace_list = [1,2,3,4,5,6,7,8];
                start_time_of_interest = start_time_of_interest_default;
                timelength_of_interest = timelength_of_interest_default;
                
            otherwise
                disp('Not ready');
                return
        end
        
        %% load data
        if isempty(filename)
            continue
        end
        desc = {[hrmin,', ',human_name,', Case#',num2str(caseid)],...
            note,[probe, ', ', imgmode,', Cc=', num2str(cc), ', MI=', num2str(mi)]};
        disp(['Loading ',data_path,filename,' ......']);
        % disp(['desc_',num2str(caseid), ' = {''',desc{1},', ',desc{2},'''};']);
        
        % read SHAPE TIC and ECG signals from csv file
        datainfo = getCSVformat([data_path,filename]);
        if ~isempty(trace_labels)
            datainfo.trace_labels = trace_labels;
        end
        
        maximum_sampling_length = 300000;
        sha_range = [datainfo.row+1,1,datainfo.row+maximum_sampling_length,datainfo.trace_indices(end)];
        sha_ts = readmatrix([data_path, filename],'Range',sha_range);
        if ~isnan(datainfo.ecg_time_index) && ~isnan(datainfo.ecg_trace_index)
            ecg_range = [datainfo.row+1,datainfo.ecg_time_index,datainfo.row+maximum_sampling_length,datainfo.ecg_trace_index];
            ecg_ts = readmatrix([data_path, filename],'Range',ecg_range);
        else
            ecg_ts = nan;
        end
        sha_ts = rmmissing(sha_ts);
        ecg_ts = rmmissing(ecg_ts);
        
        % bandstop filtering the inspiration effects
        sample_rate = 1/mean(diff(sha_ts(:,1)))
        bsFilter = designfilt('bandstopfir',...
            'FilterOrder',128,...
            'PassbandFrequency1',1.5,...
            'StopbandFrequency1',1.6,...
            'StopbandFrequency2',1.7,...
            'PassbandFrequency2',1.8,...
            'DesignMethod','ls',...
            'PassbandWeight1',1,...
            'StopbandWeight',3,...
            'PassbandWeight2',1,...
            'SampleRate',sample_rate);
        sha_ts(:,2:end) = filtfilt(bsFilter,sha_ts(:,2:end));
        
        %% define the structure for data storing
        c.code = human_name;
        c.acq_time = hrmin;
        c.note = note;
        c.caseid = caseid;
        c.filename = filename;
        c.probe = probe; %C1-6, C1-6
        c.imgmode = imgmode; %RES1, SUBAM, HRES
        c.cc = cc;
        c.mi = mi;
        c.trace = sha_ts;
        c.ecg = ecg_ts;
        c.time_index = datainfo.time_index;
        c.t1_index = datainfo.t1_index;
        c.t2_index = datainfo.t2_index;
        c.trace_indices = datainfo.trace_indices;
        c.trace_colors = datainfo.trace_colors;
        c.trace_labels = datainfo.trace_labels;
        c.visualize_trace_list = visualize_trace_list;
        c.start_time_of_interest = start_time_of_interest;
        c.timelength_of_interest = timelength_of_interest;
        c.desc = desc;
        c.ecgBeatFreq = heart_beat_freq_default; %Hz
        CASES{caseid} = c;
        
    end
    %     save(mat_file,'CASES');
else
    load(mat_file);
    disp('Load the stored case data...')
    for caseid=1:length(CASES)
        if ~isempty(CASES{caseid})
            disp(CASES{caseid}.desc);
        end
    end
end

if exist('data','var') && strcmp(data.name,human_name)
    disp('Blood Pressure Database has been loaded!');
else
    disp('Loading blood pressure database ......');
    load(bp_file);
end

%% Pre-analysis cases
for caseid = STUDY_CASE_IDS
    c = CASES{caseid};
    
    %% check the trace and ecg signals and extract for interest
    start_time_of_interest = c.start_time_of_interest;
    timelength_of_interest = c.timelength_of_interest;
    
    sha_ts = rmmissing(c.trace);
    ecg_ts = rmmissing(c.ecg);
    
    if isempty(sha_ts) || isempty(ecg_ts)
        disp(['No valid time series of subharmonic or ECG in ', c.desc{1},', ',c.desc{2}]);
        continue;
    end
    % sha_sampling_rate = (length(sha_ts(:,1))-1)/(sha_ts(end,1)-sha_ts(1,1));
    
    %% Define the time series of interest
    list_of_interest = (sha_ts(:,1)>=sha_ts(1,1)+start_time_of_interest) & ...
        (sha_ts(:,1)<(sha_ts(1,1)+start_time_of_interest+timelength_of_interest));
    sha_ts = sha_ts(list_of_interest,:);
    
    ecg_list_of_interest = (ecg_ts(:,1)>=ecg_ts(1,1)+start_time_of_interest) & ...
        (ecg_ts(:,1)<(ecg_ts(1,1)+start_time_of_interest+timelength_of_interest));
    ecg_ts = ecg_ts(ecg_list_of_interest,:);
    
    %% Identify R-waves, Pair R-wave for Beat Cycles, Normalize the time of each cycle

    [cycles, true_peak_ts] = process_sha(ecg_ts, sha_ts);

    CASES{caseid}.cycles = cycles;
    CASES{caseid}.true_peak_ts = true_peak_ts;
    CASES{caseid}.ecg_ts = ecg_ts;
    CASES{caseid}.sha_ts = sha_ts;
    
    %% Calculate the mean beat frequency
    number_of_cycles = length(cycles);
    total_cycle_time = 0;
    for cycle_ndx = 1:number_of_cycles
        total_cycle_time = total_cycle_time + cycles{cycle_ndx}.cycle_time;
    end
    CASES{caseid}.ecgBeatFreq = number_of_cycles / total_cycle_time;
end

%% Analysis & Plot the selected cases
lineStyles = {'-','-','-','-','--','--','--','--'};
lineColors = {'b','r','c','m','b','r','c','m'};
lineWidths = [1,1,1,1,1.2,1.2,1.2,1.2];

for caseid = STUDY_CASE_IDS
    
    %% load blood pressure measurement and generate the averaged signal
    c = CASES{caseid};
    
    plot_on = 1;
    [bp_time,bp_ecg,bp_bp,bp_start_time] = getBPfromDB(data,c.acq_time,c.start_time_of_interest,c.timelength_of_interest,'00:03:04');
    [bp_time_filtered,bp_bp_filtered, time_normalized, signal_normalized, signal_normalized_filtered, time_normalized_3, signal_normalized_3, signal_normalized_filtered_3] = commonCyclicSignal(bp_time,bp_bp);
    plot_on = 0;
    time_scaled4corr = 0:0.01:1;
    bp_signal_scaled4corr = interp1(time_normalized,signal_normalized_filtered,time_scaled4corr,'linear','extrap');
    
    fig=figure(1001); fig.Position=[20 844 1880 270];
    subplot(2,1,1), plot(bp_time,bp_ecg); xlabel('Time(sec)'); ylabel('ECG'); title('BP/ECG Acquisition')
    subplot(2,1,2), plot(bp_time,bp_bp); xlabel('Time(sec)'); ylabel('BP(mmHg)');
    
    %% Initiate each SHAPE case
    visualize_trace_list = c.visualize_trace_list; %[1,2,3,4,5]
    %     start_time_of_interest = c.start_time_of_interest;
    %     timelength_of_interest = c.timelength_of_interest;
    if ~isfield(c,'sha_ts') || ~isfield(c,'ecg_ts')
        disp(['Case #',num2str(caseid),' has no proper sha or ecg time series!']);
        return;
    end
    sha_ts = c.sha_ts;
    ecg_ts = c.ecg_ts;
    trace_indices = c.trace_indices;
    trace_labels = string(c.trace_labels);
    desc = c.desc;
    cycles = CASES{caseid}.cycles;
    true_peak_ts = CASES{caseid}.true_peak_ts;
    
    %% Generate scatters from all cycles and calculate the mean beat frequency
    normalized_sha_scatters = [];
    for cycle_ndx = 1:length(cycles)
        normalized_sha_scatters = cat(1,normalized_sha_scatters,cycles{cycle_ndx}.normalized_sha_ts);
    end
    normalized_sha_scatters = normalized_sha_scatters(normalized_sha_scatters(:,1)>0,:);
    
    %% Generate scatters using 3 cases (#1, #2, #3)
    % normalized_sha_scatters = [];
    % for i = 18:20
    %     for cycle_ndx = 1:length(CASES{i}.cycles)
    %         A = normalized_sha_scatters;
    %         B = CASES{i}.cycles{cycle_ndx}.normalized_sha_ts;
    %         Acols = size(A,2); Bcols = size(B,2);
    %         if Acols > Bcols && Bcols > 0
    %             B = cat(2,B, nan(size(B,1),Acols-Bcols));
    %         elseif Acols < Bcols && Acols > 0
    %             A = cat(2,A, nan(size(A,1),Bcols-Acols));
    %         end
    %         normalized_sha_scatters = cat(1,A,B);
    %     end
    % end
    % normalized_sha_scatters = normalized_sha_scatters(normalized_sha_scatters(:,1)>0,:);
    
    %% Sort scatter points
    sorted_normalized_sha_scatters = normalized_sha_scatters;
    [sorted_normalized_sha_scatters(:,1),index] = sort(normalized_sha_scatters(:,1));
    sorted_normalized_sha_scatters(:,trace_indices) = normalized_sha_scatters(index,trace_indices);
    
    %% Divide the normalized cycle into 20 sub-section to calculate the average
    interval = 1/20;
    adjusted_timeframe(:,1) = (0.0:interval:1)';
    adjusted_sha_ts = zeros(length(adjusted_timeframe)-1,trace_indices(end));
    adjusted_sha_ts(:,1) = (interval/2:interval:1)';
    for ndx=1:size(adjusted_timeframe,1)-1
        list = (sorted_normalized_sha_scatters(:,1)>adjusted_timeframe(ndx)) & ...
            (sorted_normalized_sha_scatters(:,1)<=adjusted_timeframe(ndx+1));
        for trace_ndx = trace_indices
            adjusted_sha_ts(ndx,trace_ndx) = mean(rmoutliers(sorted_normalized_sha_scatters(list,trace_ndx),'gesd'));
            % adjusted_sha_ts(ndx,trace_ndx) = mean(sorted_normalized_sha_ts(list,trace_list));
        end
    end
    
    %% Plot the raw signals of subharmonic amplitude and ECG signals
    % figure(100+caseid),
    % plotyy(sha_ts(:,1),sha_ts(:,trace_list),ecg_ts(:,1),ecg_ts(:,2));
    % legend('Yellow','L-Blue','Red','ECG');
    
    %     figure(200+caseid),
    %     plt=plotyy(sha_ts(:,1),mean(sha_ts(:,trace_list),2),ecg_ts(:,1),ecg_ts(:,2)); grid on;
    %     legend('Mean SubH','ECG');
    %     title('Mean Subharmonic Amplitude & ECG');
    %     plt(1).XLim = [sha_ts(1,1),sha_ts(end,1)];
    %     plt(2).XLim = [sha_ts(1,1),sha_ts(end,1)];
    
    %% Plot the first trace and ECG signals separately
    if caseid==1 
        traceSeq = 5;
    elseif caseid==8
        traceSeq = 1;
    else
        traceSeq = 1;
    end
    fig=figure(300+caseid);
    fig.Position=[100,100,1480,980];
    subplot(3,5,[6,7,8,9,10]);
    plot(sha_ts(:,c.time_index),sha_ts(:,c.trace_indices(traceSeq)),'-o','MarkerSize',3,'MarkerFaceColor','b');
    ylabel('Subharmonic(dB)','FontSize',14);
    xlabel('Time(s)','FontSize',14); grid on;
    title(['Subharmonic Amplitude Signal at Upstream ROI-4 on Long-Axis Image of R01'],'FontSize',14);
    xlim([sha_ts(1,1),sha_ts(end,1)]);
    subplot(3,5,[1,2,3,4,5]);
    plot(ecg_ts(:,1),ecg_ts(:,2));
    ylabel('ECG','FontSize',14);
    xlabel('Time(s)','FontSize',14); grid on;
    xlim([sha_ts(1,1),sha_ts(end,1)]);
    title('ECG Signal with R-wave Recoginized','FontSize',14)
    for cycle_ndx = 1:4
        cycle = cycles{cycle_ndx};
        fig=figure(300+caseid);
        subplot(3,5,10+cycle_ndx);
        plot(cycle.sha_ts(:,1),cycle.sha_ts(:,trace_indices(traceSeq)),'-o','MarkerSize',3,'MarkerFaceColor','b');
        title(['Cycle #', num2str(cycle_ndx)],'FontSize',14);
        xlabel('Time (s)','FontSize',14);
        ylabel('Subharmonic (dB)','FontSize',14);
        grid on
    end
    
    
    %% plot the first 5 seconds ECG and Traces
    figure(300+caseid),
    subplot(3,5,[6,7,8,9,10]); YL=ylim();
    line(ones(2,1)*true_peak_ts(:,1)',YL'*ones(1,size(true_peak_ts,1)),'LineStyle',':','Color','r','LineWidth',0.8);
    ylim(YL);
    %     xlim([start_time_of_interest,start_time_of_interest+5]);
    subplot(3,5,[1,2,3,4,5]), hold on; plot(true_peak_ts(:,1),true_peak_ts(:,2),'r*','MarkerSize',4);
    %     xlim([start_time_of_interest,start_time_of_interest+5]);
    
    %% plot the first 5 beat cycles of the first trace
    for cycle_ndx = 1:5
        cycle = cycles{cycle_ndx};
        fig=figure(400+caseid); fig.Position=[1526,42,392,1082];
        subplot(5,1,cycle_ndx);
        plot(cycle.sha_ts(:,1),cycle.sha_ts(:,trace_indices(1)),'-o','MarkerSize',3,'MarkerFaceColor','b');
        title(['Normalized Beat Cycle #', num2str(cycle_ndx)]);
        xlabel('Normalized Time [0,1]');
        ylabel('Subharmonic (dB)');
        grid on
    end
    
    %% Plot the reconstructed Subharmonic Amplitude Profile
    
    fig600_lgds = {};
    multrace_time = [];
    multrace_signal = [];
    for plotNdx=1:length(visualize_trace_list)
        if plotNdx>length(trace_indices)
            continue;
        end
        
        %SSA filtering the scatter points
        traceSeq = visualize_trace_list(plotNdx);
        time = sorted_normalized_sha_scatters(:,1);
        signal = sorted_normalized_sha_scatters(:,trace_indices(traceSeq));
        % signal = mean(sorted_normalized_sha_scatters(:,trace_indices(visualize_trace_list)),2);
        % signalFiltered = extractSignal(time, signal, 0.2);
        
        %Remove outliers
        time_outliers_removed = time;
        signal_outliers_removed = signal;
        %[signal_outliers_removed,indicators] = rmoutliers(signal,'median');
        %time_outliers_removed = time(~indicators);
        
        if ismember(traceSeq,STUDY_TRAC_IDS)
            multrace_time = cat(1,multrace_time,time_outliers_removed + traceSeq*1e-4);
            multrace_signal = cat(1,multrace_signal,signal_outliers_removed);
        end
        
        % filtering
        % signal_outliers_removed_filtered = extractSignal(time_outliers_removed, signal_outliers_removed, 0.2);
        
        % combining 3-cycle data
        combined_time = [time_outliers_removed;...
            time_outliers_removed+1;...
            time_outliers_removed+2];
        combined_signal = [signal_outliers_removed;...
            signal_outliers_removed;...
            signal_outliers_removed];
        %         combined_time = [sorted_normalized_sha_scatters(:,1);...
        %             sorted_normalized_sha_scatters(:,1)+1;...
        %             sorted_normalized_sha_scatters(:,1)+2;...
        %             sorted_normalized_sha_scatters(:,1)+3];
        %         combined_signal = [sorted_normalized_sha_scatters(:,trace_indices(traceSeq));...
        %             sorted_normalized_sha_scatters(:,trace_indices(traceSeq));...
        %             sorted_normalized_sha_scatters(:,trace_indices(traceSeq));...
        %             sorted_normalized_sha_scatters(:,trace_indices(traceSeq))];
        combined_signalFiltered = extractSignal(combined_time, combined_signal, 0.2, 1);

        resample_time = 0:0.01:3;
        resample_signal = interp1(combined_time,combined_signalFiltered,resample_time,'linear','extrap');
        
        combined_time = resample_time;
        combined_signalFiltered = resample_signal;
%         combined_signalFiltered = extractSignal(resample_time, resample_signal, 0.2, 3);
    
        second_cycle_list = combined_time>=1 & combined_time<=2;
        combined_time_2cycle = combined_time(second_cycle_list) - 1;
        combined_signalFiltered_2cycle = combined_signalFiltered(second_cycle_list);
        
        % resampling
        time_scaled4corr = combined_time_2cycle(1):0.01:combined_time_2cycle(end);
        sha_signal_scaled4corr = interp1(combined_time_2cycle,combined_signalFiltered_2cycle,time_scaled4corr,'linear','extrap');
        
        max_delay = 0; max_corr = 0; max_pval = 0;
        for ndx = 1:length(time_scaled4corr)
            x = bp_signal_scaled4corr(:);
            y = sha_signal_scaled4corr(:);
            
            if ndx~=1
                y = [y(ndx:end);y(1:ndx-1)];
            end
            
            z = [x,y]; z = rmmissing(z);
            [corr_coefs, pvalues] = corr(z);
            coef = corr_coefs(1,2); pval = pvalues(1,2);
            if coef<0 && coef<max_corr
                max_delay = ndx;
                max_corr = coef;
                max_pval = pval;
            end
        end
        [max_corr, max_pval, max_delay]
        
        %% plot the correlation between blood pressure and estimation
        fig = figure(caseid*100+plotNdx);
        x = bp_signal_scaled4corr(:);
        y = sha_signal_scaled4corr(:);
        if max_delay>1
            y = [y(max_delay:end);y(1:max_delay-1)];
        end
        min_bp = min(x); range_bp = max(x) - min(x); sensitivity = (min(y)-max(y)) / range_bp;
        z = (y - max(y)) / sensitivity + min_bp;
        plt = plot(time_scaled4corr(:),x,'b-',time_scaled4corr(:),z,'r:');
        plt(1).LineWidth = 1.2; plt(2).LineWidth = 1.6;
        legend('Measured BP','Reconstructed BP','FontSize',16);
        xlabel('Normalized Time (100%)','FontSize',16);
        ylabel('Pressure (mmHg)','FontSize',16);
        title({['corr = ',num2str(max_corr,'%2.3f'),', p-value = ',num2str(max_pval,'%2.3f')],...
            ['sensitivity= ',num2str(sensitivity,'%2.3f'),' dB/mmHg']},'FontSize',16);
        grid on;
        
        %% calculate the parameters of estiamted pressure waveform
        avgSHlst = combined_time>=1 & combined_time<2;
        avgSHAmp = mean(combined_signalFiltered(avgSHlst));
        maxSHAmp = max(combined_signalFiltered(avgSHlst));
        minSHAmp = min(combined_signalFiltered(avgSHlst));
        gapSHAmp = maxSHAmp - minSHAmp;
        
        %% plot overview of SHAPE Artery Estimation
        if (caseid==1 && traceSeq==5) || (caseid==8 && traceSeq==1)
            figure(300+caseid); subplot(3,5,15);
            plt1=plot(normalized_sha_scatters(:,1), normalized_sha_scatters(:,trace_indices(traceSeq)),'go',...
                combined_time_2cycle,combined_signalFiltered_2cycle,'r-');
            plt1(1).MarkerSize = 3; plt1(1).MarkerFaceColor='g';
            plt1(2).LineWidth = 2;
            hold on;
            plot(time_outliers_removed, signal_outliers_removed,'k.');
            title({'Projected Samples &', 'Reconstructed Waveform'},'FontSize',14);
            xlabel('Normalized Time (100%)','FontSize',14);
            ylabel('Subharmonic (dB)','FontSize',14);
            grid on;
        end
    
        fig=figure(500+caseid); fig.Position=[20 40 1880 980];
        subplot(4,length(visualize_trace_list),plotNdx),
        plt1=plot(normalized_sha_scatters(:,1), normalized_sha_scatters(:,trace_indices(traceSeq)),'k.',...
            combined_time_2cycle,combined_signalFiltered_2cycle,'r-');
        hold on;
        plot(time_outliers_removed, signal_outliers_removed,'k.');
        title(trace_labels{traceSeq});
        grid on;
        
        subplot(4,length(visualize_trace_list),length(visualize_trace_list)+plotNdx),
        plt2=plot(adjusted_sha_ts(:,1),adjusted_sha_ts(:,trace_indices(traceSeq)),'k--s',...
            combined_time_2cycle,combined_signalFiltered_2cycle,'r-'); hold on;
        %         XL = xlim(); YL = ylim(); text(XL(1),YL(1),...
        %             [num2str(avgSHAmp,'%2.1f'),'   ',num2str(maxSHAmp,'%2.1f'),'   ',num2str(minSHAmp,'%2.1f'),'   ',num2str(gapSHAmp,'%2.1f')],...
        %             'Color','r','FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
        title({['Reconstructed'],...
            ['corr = ', num2str(max_corr,'%2.3f'),', p-value = ',num2str(max_pval,'%2.3f')]});
        grid on;
        
        subplot(4,length(visualize_trace_list),length(visualize_trace_list)*2+plotNdx),
        plt3=plot(combined_time,-combined_signalFiltered,'r-'); hold on;
        YL = ylim();
        line(ones(2,1)*[1,2,3],YL'*ones(1,3),'LineStyle','-.','Color','g','LineWidth',1.0);
        ylim(YL)
        title({'3-Cycle Estimation',...
            [num2str(avgSHAmp,'%2.1f'),'   ',num2str(maxSHAmp,'%2.1f'),'   ',num2str(minSHAmp,'%2.1f'),'   ',num2str(gapSHAmp,'%2.1f')]}); grid on;
        xlabel('Normalized Time(100%)')
        
        if plotNdx==1
            ax1=subplot(4,length(visualize_trace_list),plotNdx); YL1=ylim(); hold on
            line(ax1,ones(2,1)*adjusted_timeframe',YL1'*ones(1,length(adjusted_timeframe)),...
                'LineStyle',':','Color','g','LineWidth',0.8);
            ylim(YL1);
            ylabel(ax1,{desc{1},desc{2},desc{3},'Amplitude(dB)'});
            
            ax2=subplot(4,length(visualize_trace_list),length(visualize_trace_list)+plotNdx); YL2=ylim(); hold on
            line(ax2,ones(2,1)*adjusted_timeframe',YL2'*ones(1,length(adjusted_timeframe)),...
                'LineStyle',':','Color','g','LineWidth',0.8);
            ylim(YL2);
            ylabel(ax2,{'Amplitude(dB)'});
            
            ax3=subplot(4,length(visualize_trace_list),length(visualize_trace_list)*2+plotNdx); YL3=ylim(); hold on
            line(ax3,ones(2,1)*adjusted_timeframe',YL3'*ones(1,length(adjusted_timeframe)),...
                'LineStyle',':','Color','g','LineWidth',0.8);
            ylim(YL3);
            ylabel(ax3,{'-1*Amplitude(dB)'});
            
            %         legend(plt1,{'Normalized-And-Synchronized','SSA-Filtered'});
            %         legend(plt2,{'Segmented-and-Averaged','SSA-Filtered'});
            %         legend(plt3,{'SSA-Filtered'});
        end
        
        avgBPlst = time_normalized_3>=2 & time_normalized_3<3;
        avgBPAmp = mean(signal_normalized_filtered_3(avgBPlst));
        maxBPAmp = max(signal_normalized_filtered_3(avgBPlst));
        minBPAmp = min(signal_normalized_filtered_3(avgBPlst));
        gapBPAmp = maxBPAmp - minBPAmp;
        
        subplot(4,length(visualize_trace_list),length(visualize_trace_list)*3+1),
        plot(time_normalized_3, signal_normalized_filtered_3); grid on;
        xlabel('Time(%)'); ylabel('BP Averaged (mmHg)');
        title({'3-Cycle BP',...
            [num2str(avgBPAmp,'%2.1f'),'   ',num2str(maxBPAmp,'%2.1f'),'   ',num2str(minBPAmp,'%2.1f'),'   ',num2str(gapBPAmp,'%2.1f')]});
        
        subplot(4,length(visualize_trace_list),[length(visualize_trace_list)*3+2,length(visualize_trace_list)*4]),
        plot(bp_time,bp_bp,bp_time_filtered,bp_bp_filtered); grid on;
        xlabel('Time(sec)'); ylabel('BP(mmHg)'); title('Acquired Blood Pressure')
        legend('Raw BP','Filtered BP');
        
        % plot the multiple cycles of filtered traces
        fig=figure(600+caseid); fig.Visible='off';
        plt=plot(combined_time,combined_signalFiltered,...
            'LineStyle',lineStyles{plotNdx},...
            'Color',lineColors{plotNdx},...
            'lineWidth',lineWidths(plotNdx)); hold on;
        fig600_lgds{plotNdx} = [trace_labels{traceSeq},'  ',num2str(avgSHAmp,'%2.1f'),' dB'];
        
    end
    
    %% plot all-trace-as-one analysis
    % combining 3-cycle data
    combined_time = [multrace_time;...
        multrace_time+1;...
        multrace_time+2];
    combined_signal = [multrace_signal;...
        multrace_signal;...
        multrace_signal];
    [combined_time,index] = sort(combined_time);
    combined_signal = combined_signal(index);
    
    combined_signalFiltered = extractSignal(combined_time, combined_signal, 0.2,1);
    resample_time = 0:0.01:3;
    resample_signal = interp1(combined_time,combined_signalFiltered,resample_time,'linear','extrap');
    
    combined_time = resample_time;
    combined_signalFiltered = resample_signal;
    % combined_signalFiltered = extractSignal(resample_time, resample_signal, 0.2,3);
    
    
    second_cycle_list = combined_time>=1 & combined_time<=2;
    combined_time_2cycle = combined_time(second_cycle_list) - 1;
    combined_signalFiltered_2cycle = combined_signalFiltered(second_cycle_list);
    
    % resampling
    time_scaled4corr = 0:0.01:1;
    sha_signal_scaled4corr = interp1(combined_time_2cycle,combined_signalFiltered_2cycle,time_scaled4corr,'linear','extrap');
    
    max_delay = 0; max_corr = 0; max_pval = 0;
    for ndx = 1:length(time_scaled4corr)
        x = bp_signal_scaled4corr(:);
        y = sha_signal_scaled4corr(:);
        
        if ndx~=1
            y = [y(ndx:end);y(1:ndx-1)];
        end
        
        z = [x,y]; z = rmmissing(z);
        [corr_coefs, pvalues] = corr(z);
        coef = corr_coefs(1,2); pval = pvalues(1,2);
        if coef<0 && coef<max_corr
            max_delay = ndx;
            max_corr = coef;
            max_pval = pval;
        end
    end
    [max_corr, max_pval, max_delay]
    
    %% plot the correlation between blood pressure and estimation
    fig = figure(caseid*1000); fig.Position = [800,120,560,420];
    x = bp_signal_scaled4corr(:);
    y = sha_signal_scaled4corr(:);
    if max_delay>1
        y = [y(max_delay:end);y(1:max_delay-1)];
    end
    min_bp = min(x); range_bp = max(x) - min(x); sensitivity = (min(y)-max(y)) / range_bp;
    z = (y - max(y)) / sensitivity + min_bp;
    plt = plot(time_scaled4corr(:),x,'b-',time_scaled4corr(:),z,'r:');
    plt(1).LineWidth = 1; plt(2).LineWidth = 1.4;
    legend('Measured BP','Reconstructed BP');
    xlabel('Normalized Time (100%)');
    ylabel('Pressure (mmHg)');
    title({['Comparison of Reconstructed and Measured BP'],...
        ['corr = ',num2str(max_corr,'%2.3f'),', p-value = ',num2str(max_pval,'%2.3f')],...
        ['sensitivity= ',num2str(sensitivity,'%2.3f'),' dB/mmHg']});
    grid on;
    
    %% calculate the parameters of estiamted pressure waveform
    avgSHlst = combined_time>=1 & combined_time<2;
    avgSHAmp = mean(combined_signalFiltered(avgSHlst));
    maxSHAmp = max(combined_signalFiltered(avgSHlst));
    minSHAmp = min(combined_signalFiltered(avgSHlst));
    gapSHAmp = maxSHAmp - minSHAmp;
    
    %% plot overview of SHAPE Artery Estimation
    fig=figure(5000+caseid); fig.Position=[100 200 480 800];
    subplot(4,1,1),
    plt1=plot(multrace_time, multrace_signal,'k.',...
        combined_time_2cycle,combined_signalFiltered_2cycle,'r-');hold on;grid on;
    title('Scatters of All Traces');
    line(ax1,ones(2,1)*adjusted_timeframe',YL1'*ones(1,length(adjusted_timeframe)),...
        'LineStyle',':','Color','g','LineWidth',0.8);
    ylim(YL1);
    ylabel(ax1,{desc{1},desc{2},desc{3},'Amplitude(dB)'});
    
    subplot(4,1,2),
    plt2=plot(combined_time_2cycle,combined_signalFiltered_2cycle,'r-'); hold on;grid on;
    title({['Reconstructed'],...
        ['corr = ', num2str(max_corr,'%2.3f'),', p-value = ',num2str(max_pval,'%2.3f')]});
    YL2=ylim();
    line(ax2,ones(2,1)*adjusted_timeframe',YL2'*ones(1,length(adjusted_timeframe)),...
        'LineStyle',':','Color','g','LineWidth',0.8);
    ylim(YL2);
    ylabel(ax2,{'Amplitude(dB)'});
    
    subplot(4,1,3),
    plt3=plot(combined_time,-combined_signalFiltered,'r-'); hold on;
    YL3 = ylim();
    line(ones(2,1)*[1,2,3],YL3'*ones(1,3),'LineStyle','-.','Color','g','LineWidth',1.0);
    title({'3-Cycle Estimation',...
        [num2str(avgSHAmp,'%2.1f'),'   ',num2str(maxSHAmp,'%2.1f'),'   ',num2str(minSHAmp,'%2.1f'),'   ',num2str(gapSHAmp,'%2.1f')]}); grid on;
    xlabel('Normalized Time(100%)')
    ylabel(ax3,{'-1*Amplitude(dB)'});
    ylim(YL3);
    
    avgBPlst = time_normalized_3>=2 & time_normalized_3<3;
    avgBPAmp = mean(signal_normalized_filtered_3(avgBPlst));
    maxBPAmp = max(signal_normalized_filtered_3(avgBPlst));
    minBPAmp = min(signal_normalized_filtered_3(avgBPlst));
    gapBPAmp = maxBPAmp - minBPAmp;
    
    subplot(4,1,4),
    plot(time_normalized_3, signal_normalized_filtered_3); grid on;
    xlabel('Time(%)'); ylabel('BP Averaged (mmHg)');
    title({'3-Cycle BP',...
        [num2str(avgBPAmp,'%2.1f'),'   ',num2str(maxBPAmp,'%2.1f'),'   ',num2str(minBPAmp,'%2.1f'),'   ',num2str(gapBPAmp,'%2.1f')]});
    
    
    % plot the multiple cycles of filtered traces
    fig=figure(600+caseid); fig.Visible='off';
    plt=plot(combined_time,combined_signalFiltered,...
        'LineStyle',lineStyles{plotNdx},...
        'Color',lineColors{plotNdx},...
        'lineWidth',lineWidths(plotNdx)); hold on;
    fig600_lgds{plotNdx} = [trace_labels{traceSeq},'  ',num2str(avgSHAmp,'%2.1f'),' dB'];
    
end

