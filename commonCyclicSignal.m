function [time_filtered,signal_filtered, time_normalized, signal_normalized, signal_normalized_filtered, time_normalized_3, signal_normalized_3, signal_normalized_filtered_3] = commonCyclicSignal(tm1,bp1,fs)
% identify and pair the peak and valley of cyclic signal,
% synchronize the cycles and plot the scatters of each cycle
% extract the mean signal from teh synchronized scatters
if nargin<3
    fs = 1000;
end

global plot_on

downsample_rate= round(fs / 250);
tm = downsample(tm1, downsample_rate); %100Hz; 降采样到100Hz
bp = downsample(bp1, downsample_rate); %100Hz; 降采样到100Hz
fs = fs / downsample_rate;

%% 频谱分析，提取心率
N = length(bp);
y = fft(bp-mean(bp),N);
f = (0:N-1)'*fs/N;
fLen = round(N/2);
frqs = f(1:fLen);
mags = abs(y(1:fLen));
if plot_on
    figure(); 
    subplot(2,1,1), plot(tm1,bp1);
    xlabel('Time(s)');
    ylabel('Raw Blood Pressure (mmHg)');
    subplot(2,1,2), plot(frqs(2:end),mags(2:end));
    xlabel('Freq(Hz)'); ylabel('Mag.')
end

mags(frqs<0.5 | frqs>6) = 0;
[MaxMag,MaxId] = max(mags);
beat_freq = frqs(MaxId(1));

%% 信号滤波，剔除f/10~f/2Hz范围的信号
L=240; Np=20; frqsRemoved = [beat_freq/20,beat_freq/2]; prcpRemoved = []; prcpRemained = [];
[signal_extracted] = ssa(bp,L,Np);

if plot_on
    fig=figure(); fig.Position=[962 42 958 1078];
end
plotNum = Np;
for i=1:plotNum
    x = signal_extracted(:,i);
    %fs = 50; %50Hz
    N = length(x);
    y = fft(x,N);
    f = (0:N-1)'*fs/N;
    fLen = round(N/2);
    frqs = f(1:fLen);
    mags = abs(y(1:fLen));
    
    maxMag = max(mags);
    frqLst = find((frqs>frqsRemoved(1)) & (frqs<frqsRemoved(2)));
    maxSct = max(mags(frqLst));
    if maxSct > maxMag * 0.2
        prcpRemoved = [prcpRemoved, i];
    else
        prcpRemained = [prcpRemained, i];
    end
    
    if plot_on
        figure(fig);
        subplot(plotNum,2,2*i-1), plot(tm,x,'b-'); grid on; ylabel(['P',num2str(i)]);
        subplot(plotNum,2,2*i), plot(frqs,mags,'b-'); grid on; ylabel(['M',num2str(i)]);
    end
end
if plot_on
    disp('Beat Frequency:');
    disp(beat_freq);
    disp('Removed Principal Components : ');
    disp(prcpRemoved);
    disp('Retained Principal Components : ');
    disp(prcpRemained);
end

principles = {prcpRemoved, prcpRemained};
Xall = sum(signal_extracted,2);
bp_removed  = sum(signal_extracted(:,principles{1}),2);
% if principles{2}(1)==1
%     bp_remained  = sum(signalExtracted(:,principles{2}(2:end)),2) + mean(signalExtracted(:,principles{2}(1)));
% else
    bp_remained  = sum(signal_extracted(:,principles{2}),2);
% end

% bsFilt = designfilt('bandstopfir','FilterOrder',400,...
%     'PassbandFrequency1',0.02,'PassbandFrequency2',1.0,...
%     'StopbandFrequency1',0.15,'StopbandFrequency2',0.5,...
%     'DesignMethod','ls',...
%     'SampleRate',50);
% fvtool(bsFilt)
% bp_filtered = filtfilt(bsFilt,Xall);

if plot_on
    fig=figure(); fig.Position=[7 667 1899 420];
    subplot(1,2,1), plot(tm,bp,'b-',tm,Xall,'r:'); grid on
    subplot(1,2,2), plot(tm,bp_removed,'b-',tm,bp_remained,'r-'); grid on
end

%% segment the bp_remained to get the average bp profile
signal_filtered = bp_remained;
time_filtered = tm;
timeInterval = 2/fs; %连续3个采样点

gradient = nan(size(signal_filtered));
for ndx = 1:length(time_filtered)-1
    currentTime = time_filtered(ndx);
%     rabList = (timesq>(currentTime-timeInterval)) & (timesq<(currentTime+timeInterval)); %取前后相邻共5个采样点计算梯度
%     st = timesq(rabList); sx = signal(rabList);
%     if length(st)>=3
%         [poly_coef,~] = polyfit(st, sx, 2);
%         gradient(ndx) = 2*poly_coef(2)*currentTime+poly_coef(1); %the first coefficient is the gradient
%     end
    gradient(ndx) = (signal_filtered(ndx+1)-signal_filtered(ndx)) / (time_filtered(ndx+1) - time_filtered(ndx));
end

max_gradient = max(gradient,[],'all','omitnan');
min_gradient = min(gradient,[],'all','omitnan');
max_signal = max(signal_filtered,[],'all','omitnan');
min_signal = min(signal_filtered,[],'all','omitnan');

half_winL = round(fs / beat_freq / 2);
valleyT = nan(500,1);
peakT = nan(500,1);
valleyT_ndx = 0; peakT_ndx = 0;
for ndx = 1:length(time_filtered)-1 %找局部峰值或谷值点
    local_peak = 0; local_valley = 0;
    id_L = max(1,ndx-half_winL);
    id_R = min(ndx+half_winL, length(time_filtered));
    if ndx<2 || isnan(gradient(ndx-1)) || isnan(gradient(ndx)) || isnan(gradient(ndx+1))
        continue
    end
    
    if abs(signal_filtered(ndx)-max(signal_filtered(id_L:id_R)))<1e-6
        local_peak = 1;
    else
        local_peak = 0;
    end
    if abs(signal_filtered(ndx)-min(signal_filtered(id_L:id_R)))<1e-6
        local_valley = 1;
    else
        local_valley = 0;
    end
    
    if local_valley==1
        if abs(gradient(ndx))<=(max_gradient-min_gradient)*1e-6
            valleyT_ndx = valleyT_ndx + 1;
            valleyT(valleyT_ndx) = time_filtered(ndx);
        elseif (gradient(ndx-1) * gradient(ndx+1) <0)
            valleyT_ndx = valleyT_ndx + 1;
            t1 = time_filtered(ndx-1); t2 = time_filtered(ndx+1);
            x1 = gradient(ndx-1); x2 = gradient(ndx+1);
            valleyT(valleyT_ndx) = (x1*t2-x2*t1)/(x1-x2);
        end
    elseif local_peak==1
        if abs(gradient(ndx))<=(max_gradient-min_gradient)*1e-6
            peakT_ndx = peakT_ndx + 1;
            peakT(peakT_ndx) = time_filtered(ndx);
        elseif (gradient(ndx-1) * gradient(ndx+1) <0)
            peakT_ndx = peakT_ndx + 1;
            t1 = time_filtered(ndx-1); t2 = time_filtered(ndx+1);
            x1 = gradient(ndx-1); x2 = gradient(ndx+1);
            peakT(peakT_ndx) = (x1*t2-x2*t1)/(x1-x2);
        end
    else
        continue
    end
end
peakT = rmmissing(peakT);
valleyT = rmmissing(valleyT);

% pairing the peak-valley
valleyT_ndx = 1; peakT_ndx = 1;
valley_peak = nan(500,2); ndx = 0;
while valleyT_ndx<length(valleyT) && peakT_ndx<length(peakT)
    if valleyT(valleyT_ndx) < peakT(peakT_ndx)
        if (valleyT(valleyT_ndx+1) > peakT(peakT_ndx)) && ((peakT(peakT_ndx)-valleyT(valleyT_ndx))<0.7/beat_freq)
            ndx = ndx + 1;
            valley_peak(ndx,:) = [valleyT(valleyT_ndx),peakT(peakT_ndx)];
            valleyT_ndx = valleyT_ndx + 1;
            peakT_ndx = peakT_ndx + 1;
        else
            valleyT_ndx = valleyT_ndx + 1;
        end
    elseif valleyT(valleyT_ndx) > peakT(peakT_ndx)
        peakT_ndx = peakT_ndx + 1;
    else
        disp(['Incorrect peak-valley configuration at valley(',num2str(valleyT_ndx),')']);
    end
end
valley_peak = rmmissing(valley_peak);

% validate the pairs
pairing = nan(500,2); ndx = 0;
valley_ndx = 1; valley = valley_peak(:,1);
while valley_ndx<length(valley)
    gap = valley(valley_ndx+1) - valley(valley_ndx);
    if gap>0.8/beat_freq && gap<1.2/beat_freq
        ndx = ndx + 1;
        pairing(ndx,:) = [valley(valley_ndx), valley(valley_ndx+1)];
    end
    valley_ndx = valley_ndx + 1;
end
pairing = rmmissing(pairing);
cyclicValleyTime = unique([pairing(:,1);pairing(:,2)]);
% disp(cyclicValleyTime');

% normalized the time of cycles
time_sync = nan(size(time_filtered));
signal_sync = nan(size(signal_filtered));
for ndx=1:size(pairing,1)
    rabList = time_filtered>pairing(ndx,1) & time_filtered<pairing(ndx,2);
    time_sync(rabList) = (time_filtered(rabList) - pairing(ndx,1))./(pairing(ndx,2)-pairing(ndx,1));
    signal_sync(rabList) = signal_filtered(rabList);
end
% time_sync = nan(size(time_filtered));
% signal_sync = nan(size(signal_extracted));
% for ndx=1:size(pairing,1)
%     rabList = time_filtered>pairing(ndx,1) & time_filtered<pairing(ndx,2);
%     time_sync(rabList) = (time_filtered(rabList) - pairing(ndx,1))./(pairing(ndx,2)-pairing(ndx,1));
%     signal_sync(rabList) = signal_extracted(rabList);
% end

% re-align the time and filter the signal
[sorted_time_sync,indices] = sort(time_sync);
sorted_signal_sync = signal_sync(indices);
time_normalized = rmmissing(sorted_time_sync);
signal_normalized = rmmissing(sorted_signal_sync);
if ~isempty(signal_normalized)
    signal_normalized_filtered = extractSignal(time_normalized,signal_normalized,0.05);
else
    signal_normalized_filtered = signal_normalized;
end

time_normalized_3 = [time_normalized; time_normalized+1; time_normalized+2];
signal_normalized_3 = [signal_normalized; signal_normalized; signal_normalized];
if ~isempty(signal_normalized_3)
    signal_normalized_filtered_3 = extractSignal(time_normalized_3,signal_normalized_3,0.05);
else
    signal_normalized_filtered_3 = signal_normalized_3;
end

if plot_on
    fig=figure(); fig.Position = [10 80 1480 980];
    
    subplot(3,3,[1,2,3]), plot(tm1,bp1,'b-',tm,bp,'r:'); grid on; hold on;YL = ylim();
%     subplot(3,3,[1,2,3]), plot(tm1,bp1,'b-'); grid on; hold on;YL = ylim();
    %     disp([cyclicValleyTime';cyclicValleyTime'])
    %     disp(YL'*ones(1,length(cyclicValleyTime)))
    %     line([cyclicValleyTime';cyclicValleyTime'],YL'*ones(1,length(cyclicValleyTime)),'Color','g','LineStyle','-','LineWidth',0.75);
    xlim([tm(1), tm(end)])
    xlabel('Time (s)');
    legend('Raw BP','Down-Sampled BP')
    ylabel('Blood Pressure (mmHg)');
    title('Raw Signal of Blood Pressure');
    
    subplot(3,3,[4,5,6]), plt=plot(time_filtered,signal_filtered,'b-'); grid on; hold on; YL = ylim();
    line([cyclicValleyTime';cyclicValleyTime'],YL'*ones(1,length(cyclicValleyTime)),'Color','r','LineStyle',':','LineWidth',0.75);
    xlim([tm(1), tm(end)])
    xlabel('Time (s)');
    ylabel('Blood Pressure (mmHg)');
    title('Bandpass Filtered Signal of Blood Pressure ([f/20,f/2])')
    
    subplot(3,3,7)
    plt=plot(time_normalized,signal_normalized,'b.',time_normalized,signal_normalized_filtered,'r-'); grid on;
    plt(2).LineWidth = 2;
    xlabel('Normalized Time (100%)'); ylabel('Blood Pressure (mmHg)'); title('Cyclic BP Signals');
    
    subplot(3,3,8)
    plt=plot(time_normalized,signal_normalized_filtered,'r-'); grid on;
    xlabel('Normalized Time (100%)'); ylabel('Blood Pressure (mmHg)'); title('Reconstructed One-Cycle BP Waveform');
    
    subplot(3,3,9)
    plt=plot([time_normalized_3],[signal_normalized_filtered_3],'r-'); grid on;
    xlabel('Normalized Time (100%)'); ylabel('Blood Pressure (mmHg)'); title('Repetitive BP Waveform');
    
    fig = figure(); fig.Position = [60,400,1640,600];
    subplot(2,5,[1,2,3]), plot(tm1,bp1,'b-'); grid on; hold on;YL = ylim();
    %     disp([cyclicValleyTime';cyclicValleyTime'])
    %     disp(YL'*ones(1,length(cyclicValleyTime)))
    %     line([cyclicValleyTime';cyclicValleyTime'],YL'*ones(1,length(cyclicValleyTime)),'Color','g','LineStyle','-','LineWidth',0.75);
    xlim([tm(1), tm(end)])
    xlabel('Time (s)','FontSize',14);
    legend('Raw BP','FontSize',14)
    ylabel('Blood Pressure (mmHg)','FontSize',14);
    title('Raw Signal of Blood Pressure','FontSize',14);
    
    subplot(2,5,[6,7,8]), plt=plot(time_filtered,signal_filtered,'b-'); grid on; hold on; YL = ylim();
    line([cyclicValleyTime';cyclicValleyTime'],YL'*ones(1,length(cyclicValleyTime)),'Color','r','LineStyle',':','LineWidth',0.75);
    xlim([tm(1), tm(end)])
    xlabel('Time (s)','FontSize',14);
    ylabel('Blood Pressure (mmHg)','FontSize',14);
    title('Bandpass Filtered ([f/20,f/2]) BP Signal with Cardiac Cycles Recognized','FontSize',14)
    
    subplot(2,5,[4,5,9,10])
    plt=plot(time_normalized,signal_normalized,'b.',time_normalized,signal_normalized_filtered,'r-'); grid on;
    plt(2).LineWidth = 2;
    xlabel('Normalized Time (100%)','FontSize',14); 
    ylabel('Blood Pressure (mmHg)','FontSize',14); 
    title('BP Signal Reconstruction','FontSize',14);
    legend({'Projected BP Samples','Reconstructed BP Waveform'},'FontSize',14);
    
    
    %     fig=figure(6); fig.Position = [664 154 1120 840];
    %     subplot(2,1,1), plot(timesq,signal,'b-',timesq,signal,'b.'); xlabel('Time(s)'); ylabel('Pressure(mmHg)'); grid on;
    %     hold on; YL = ylim();
    %     line([pairing(:,1),pairing(:,1)]',(ones(length(pairing(:,1)),1)*YL)','LineStyle','-','Color','g');
    %     % line([valley_peak(:,2),valley_peak(:,2)]',(ones(length(valley_peak(:,2)),1)*YL)','LineStyle','-','Color','r');
    %     subplot(2,1,2), plot(timesq,gradient,'b-s'); xlabel('Time(s)'); ylabel('Gradient(mmHg/s)'); grid on;
    
    %     fig=figure(7); fig.Position = [160 160 1680 420];
    %     subplot(1,3,1)
    %     plt=plot(time_processed,signal_processed,'b.',time_processed,signal_filtered,'r-'); grid on;
    %     plt(2).LineWidth = 2;
    %     xlabel('Normalized Time(100%)'); ylabel('Blood Pressure(mmHg)');
    %
    %     subplot(1,3,2)
    %     plt=plot(time_processed,signal_filtered,'r-'); grid on;
    %     xlabel('Normalized Time(100%)'); ylabel('Averaged Blood Pressure(mmHg)');
    %
    %     subplot(1,3,3)
    %     plt=plot([time_processed_3],[signal_filtered_3],'r-'); grid on;
    %     xlabel('Normalized Time(100%)'); ylabel('Averaged Blood Pressure(mmHg)');
end
