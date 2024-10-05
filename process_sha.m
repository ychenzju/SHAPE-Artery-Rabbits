function [cycles, true_peak_ts] = process_sha(ecg_ts, sha_ts)

global plot_on

%% R-wave identification
[true_peak_ts] = procee_ecg(ecg_ts,0.4);

estimated_ecg_cycle_time = (true_peak_ts(end,1) - true_peak_ts(1,1)) / (size(true_peak_ts,1) - 1);
ecg_beat_frequency = 1 / estimated_ecg_cycle_time;

% % frequency-domain approach to estimate the beat frequency
% % it looks not like a robust method.
% fs = 500; T= ecg_ts(:,1); X = ecg_ts(:,2) - mean(ecg_ts(:,2));
% N = length(X);
% y = fft(X,N);
% f = (0:N-1)'*fs/N;
% fLen = round(N/2);
% frqs = f(1:fLen);
% mags = abs(y(1:fLen));
% if plot_on
%     figure(); 
%     subplot(2,1,1), plot(T,X);
%     xlabel('Time(s)');
%     ylabel('ECG');
%     subplot(2,1,2), plot(frqs(1:end),mags(1:end));
%     xlabel('Freq(Hz)'); ylabel('Mag.')
% end
% mags(frqs<0.5 | frqs>6) = 0;
% [~,MaxId] = max(mags);
% beat_freq = frqs(MaxId(1));
% if plot_on
%     disp('Time and Freq-Domain Estimated: ');
%     disp([ecg_beat_frequency, beat_freq]);
% end

%% normalize the beat cycles to formulate one beat-cycle subharmonic amplitude signal
len_true_peaks = size(true_peak_ts,1);
cycles = {}; cycle_ndx = 0;
for ndx=1:len_true_peaks-1
    cycle_start_time = true_peak_ts(ndx,1);
    cycle_stop_time = true_peak_ts(ndx+1,1);
    if (cycle_stop_time-cycle_start_time)>1.25/ecg_beat_frequency
        continue;
    end
    if (cycle_stop_time-cycle_start_time)<0.75/ecg_beat_frequency
        continue;
    end
    cycle_ndx= cycle_ndx + 1;
    cycle_time = cycle_stop_time - cycle_start_time;
    
    sha_list = (sha_ts(:,1)>cycle_start_time) & (sha_ts(:,1)<=cycle_stop_time);
    ecg_list = (ecg_ts(:,1)>cycle_start_time) & (ecg_ts(:,1)<=cycle_stop_time);
    normalized_ecg_ts = ecg_ts(ecg_list,:);
    normalized_ecg_ts(:,1) = (normalized_ecg_ts(:,1) - cycle_start_time) / cycle_time;
    normalized_sha_ts = sha_ts(sha_list,:);
    normalized_sha_ts(:,1) = (normalized_sha_ts(:,1) - cycle_start_time) / cycle_time;

    cycle.ecg_ts = ecg_ts(ecg_list,:);
    cycle.sha_ts = sha_ts(sha_list,:);
    cycle.normalized_sha_ts = normalized_sha_ts; 
    cycle.normalized_ecg_ts = normalized_ecg_ts;
    cycle.cycle_time = cycle_time;
    cycles{cycle_ndx} = cycle;
end

% %% generate scatters from all cycles and calculate the mean beat frequency
% normalized_sha_scatters = [];
% number_of_cycles = length(cycles);
% total_cycle_time = 0;
% for cycle_ndx = 1:number_of_cycles
%     total_cycle_time = total_cycle_time + cycles{cycle_ndx}.cycle_time;
%     normalized_sha_scatters = cat(1,normalized_sha_scatters,cycles{cycle_ndx}.normalized_sha_ts);
% end
% normalized_sha_scatters = normalized_sha_scatters(normalized_sha_scatters(:,1)>0,:);
% ecg_beat_frequency = number_of_cycles / total_cycle_time;

end

