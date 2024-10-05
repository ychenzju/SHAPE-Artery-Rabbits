function [true_peak_ts] = procee_ecg(ecg_ts,gradient_threshold_ratio)

if nargin<2
    gradient_threshold_ratio = 0.5;
end

% load debug_process_ecg_data

global plot_on

%% parameter configuration
positive_gradient_threshold = 4000; %to recognize R wave
negative_gradient_threshold = -6000;
maximum_beat_cycles = 500; %to store the information of beat cycles
interval_threshold_btw_ascending_and_descending = 30; %to pair the ascending and descending points
peak_fitting_window = 3; % the data length for peak curve fitting = 2 * peak_fitting_window + 1
peak_fitting_order = 3; % the order of polynomial curve fitting for peak point

%% calculating the ecg gradients, for identifying R waves
ecg_gradients = ecg_ts;
ecg_gradients_abs = ecg_ts;
for ndx = 1:size(ecg_gradients,1)
    if ndx==1
        ecg_gradients(ndx,2) = (ecg_ts(ndx+1,2) - ecg_ts(ndx,2)) / (ecg_ts(ndx+1,1) - ecg_ts(ndx,1));
    elseif ndx==size(ecg_gradients,1)
        ecg_gradients(ndx,2) = (ecg_ts(ndx,2) - ecg_ts(ndx-1,2)) / (ecg_ts(ndx,1) - ecg_ts(ndx-1,1));
    else
        ecg_gradients(ndx,2) = (ecg_ts(ndx+1,2) - ecg_ts(ndx-1,2)) / (ecg_ts(ndx+1,1) - ecg_ts(ndx-1,1));
    end
    ecg_gradients_abs(ndx,2) = abs(ecg_gradients(ndx,2));
end


%% calculating the proper gradient thresholds for R wave identification
maxSetsCount = 0; Knum = 0;
while maxSetsCount<5
    Knum = Knum + 3;
    maxK_gradient = maxk(ecg_gradients(:,2),Knum);
    maxSetsCount = count_sets(ecg_gradients(:,2)>maxK_gradient(end));
end
minSetsCount = 0; Knum = 1;
while minSetsCount<5
    Knum = Knum + 3;
    minK_gradient = mink(ecg_gradients(:,2),Knum);
    minSetsCount = count_sets(ecg_gradients(:,2)<minK_gradient(end));
end

max_gradient = maxK_gradient(end);
min_gradient = minK_gradient(end);
positive_gradient_threshold = gradient_threshold_ratio * max_gradient(1);
negative_gradient_threshold = gradient_threshold_ratio * min_gradient(1);

%% ideintifying the steepest ascending and descending points of R waves
ecg_gradients_gt_p10000 = (ecg_gradients(:,2)>positive_gradient_threshold);
ecg_gradients_lt_n10000 = (ecg_gradients(:,2)<negative_gradient_threshold);

sets_p10000 = count_sets(ecg_gradients_gt_p10000);
sets_n10000 = count_sets(ecg_gradients_lt_n10000);

threshold_not_found = 1;
while threshold_not_found
   
    if plot_on
        disp('ECG R Identification -- Thresholds for positive/negative gradients');
        disp([sets_p10000, sets_n10000, positive_gradient_threshold, negative_gradient_threshold]);
    end
    
    if sets_n10000>2.0*sets_p10000-3
        % 异常Q波，导致R波前缘存在大的负梯度
        threshold_not_found = 0;
    elseif sets_p10000>2.0*sets_n10000-3
        % 异常R波，导致R波后缘存在大的正梯度
        threshold_not_found = 0;
    elseif sets_n10000>=sets_p10000+3
        % 异常Q波，导致R波前缘存在大的负梯度
        negative_gradient_threshold = negative_gradient_threshold * 0.4;
    elseif (sets_p10000>=sets_n10000+3)
        negative_gradient_threshold = negative_gradient_threshold * 0.75;
    else
        threshold_not_found = 0;
    end

    if threshold_not_found==1
        ecg_gradients_gt_p10000 = (ecg_gradients(:,2)>positive_gradient_threshold);
        ecg_gradients_lt_n10000 = (ecg_gradients(:,2)<negative_gradient_threshold);
        
        sets_p10000 = count_sets(ecg_gradients_gt_p10000);
        sets_n10000 = count_sets(ecg_gradients_lt_n10000);
        
        if negative_gradient_threshold>-1500 || positive_gradient_threshold<3000
            threshold_not_found = 0;
        end
    end

end

if plot_on==1
    fig1=figure;
    fig1.Position=[620 80 1200 1000];
    subplot(3,1,1), plot(ecg_ts(:,1),ecg_ts(:,2),'b-'); ylabel('ECG'); xlabel('Time(s)');
    subplot(3,1,2), plot(ecg_gradients(:,1),ecg_gradients(:,2),'b-'); ylabel('ECG Gradient'); xlabel('Time(s)'); 
    hold on
    XL = xlim(); YL = ylim();
    line([XL(1),XL(2)],[positive_gradient_threshold,positive_gradient_threshold],'Color','red','LineStyle','--');
    line([XL(1),XL(2)],[negative_gradient_threshold,negative_gradient_threshold],'Color','green','LineStyle','--');
    xlim([ecg_ts(1,1),ecg_ts(end,1)]);
    hold off
    subplot(3,1,3), plot(ecg_gradients_abs(:,1),ecg_gradients_abs(:,2),'b-'); ylabel('ECG Gradient Abs'); xlabel('Time(s)');
    hold on
    XL = xlim(); YL = ylim();
    line([XL(1),XL(2)],[positive_gradient_threshold,positive_gradient_threshold],'Color','red','LineStyle','--');
    line([XL(1),XL(2)],[-negative_gradient_threshold,-negative_gradient_threshold],'Color','green','LineStyle','--');
    xlim([ecg_ts(1,1),ecg_ts(end,1)]);
    hold off
end

maximum_indices = zeros(maximum_beat_cycles,1);
minimum_indices = zeros(maximum_beat_cycles,1);
maximum_id = 0;
minimum_id = 0;

%% 遍历整个梯度序列，找出大于正值阈值的所有局部最大的梯度
ndx = 1;
while ndx<size(ecg_gradients,1)
    start_ndx = 0;
    stop_ndx = 0;
    if ecg_gradients_gt_p10000(ndx)==1
        %search for maximal positive gradient
        start_ndx = ndx;
        ndx = ndx + 1;
        while ecg_gradients_gt_p10000(ndx)==1 && ndx<size(ecg_gradients,1)
            ndx = ndx + 1;
        end
        stop_ndx = ndx-1;
    else
        ndx = ndx + 1;
    end
    
    if start_ndx>0 && stop_ndx>0 && start_ndx<=stop_ndx
        indices = start_ndx:stop_ndx;
        [~, id] = max(ecg_gradients(indices,2));
        maximum_ndx = start_ndx - 1 + id;
        maximum_id = maximum_id + 1;
        maximum_indices(maximum_id) = maximum_ndx;
    end
end

%% 遍历整个梯度序列，找出所有小于负值阈值的所有局部最小梯度
ndx = 1;
while ndx<size(ecg_gradients,1)
    start_ndx = 0;
    stop_ndx = 0;
    if ecg_gradients_lt_n10000(ndx)==1
        %search for minimal negative gradient
        start_ndx = ndx;
        ndx = ndx + 1;
        while ecg_gradients_lt_n10000(ndx)==1 && ndx<size(ecg_gradients,1)
            ndx = ndx + 1;
        end
        stop_ndx = ndx - 1;
    else
        ndx = ndx + 1;
    end
    
    if start_ndx>0 && stop_ndx>0 && start_ndx<=stop_ndx
        indices = start_ndx:stop_ndx;
        [~, id] = min(ecg_gradients(indices,2));
        minimum_ndx = start_ndx - 1 + id;
        minimum_id = minimum_id + 1;
        minimum_indices(minimum_id) = minimum_ndx;
    end
end
maximum_indices = maximum_indices(maximum_indices>0);
minimum_indices = minimum_indices(minimum_indices>0);

if plot_on==1
    figure(fig1);
    subplot(3,1,1), 
    hold on;
    plot(ecg_ts(maximum_indices,1),ecg_ts(maximum_indices,2),'rs');
    plot(ecg_ts(minimum_indices,1),ecg_ts(minimum_indices,2),'gs');
    xlim([ecg_ts(1,1),ecg_ts(end,1)]);
    hold off;
    
%     fig=figure(9002);fig.Position=[40 60 1078 1002];
%     subplot(3,1,1), plot(ecg_ts(:,1),ecg_ts(:,2),'b-'); ylabel('ECG');
%     hold on;
%     plot(ecg_ts(maximum_indices,1),ecg_ts(maximum_indices,2),'rs');
%     plot(ecg_ts(minimum_indices,1),ecg_ts(minimum_indices,2),'gs');
%     hold off;
%     subplot(3,1,2), plot(ecg_gradients(:,1),ecg_gradients(:,2),'b-'); ylabel('ECG Gradient')
%     subplot(3,1,3), plot(ecg_gradients(:,1),ecg_gradients_gt_p10000,ecg_gradients(:,1),ecg_gradients_lt_n10000); ylabel('frames');
end

%% pairing the ascending and descending points
ndx_maximum = 1; len_maximum = length(maximum_indices);
ndx_minimum = 1; len_minimum = length(minimum_indices);
pairing_indices = zeros(min(len_maximum,len_minimum),2);
ndx = 1;
while ndx_maximum<=len_maximum
    while ndx_minimum<=len_minimum
        ndx_gap = minimum_indices(ndx_minimum) - maximum_indices(ndx_maximum);
        if ndx_gap < 0
            ndx_minimum = ndx_minimum + 1;
        elseif ndx_gap < interval_threshold_btw_ascending_and_descending
            pairing_indices(ndx,:) = [maximum_indices(ndx_maximum), minimum_indices(ndx_minimum)];
            ndx_maximum = ndx_maximum + 1;
            ndx_minimum = ndx_minimum + 1;
            ndx = ndx + 1;
        else
            ndx_maximum = ndx_maximum + 1;
        end
        if ndx_maximum>len_maximum
            break;
        end
    end
    if ndx_minimum>len_minimum || ndx_maximum>len_maximum
        break
    end
end
pairing_indices = pairing_indices(pairing_indices(:,1)>0,:);

if plot_on==1
    fig2=figure; fig2.Position=[40 680 1750 420];
    title('Pairing and Peak Identification')
    plot(ecg_ts(:,1),ecg_ts(:,2),'b-'); hold on; ylabel('ECG'); xlabel('Time(s)');
    plot(ecg_ts(pairing_indices(:,1),1),ecg_ts(pairing_indices(:,1),2),'r.');
    plot(ecg_ts(pairing_indices(:,2),1),ecg_ts(pairing_indices(:,2),2),'g.');
    xlim([ecg_ts(1,1),ecg_ts(end,1)]);
end

%% searching for the peak points, fitting the peak curves and identify the true peak points
len_peaks = size(pairing_indices,1);
true_peak_ts = zeros(len_peaks,2);
virtual_peak_ts = zeros(len_peaks,2);
for ndx = 1:len_peaks
    [~,n] = max(ecg_ts(pairing_indices(ndx,1):pairing_indices(ndx,2),2));
    virtual_peak_ndx = pairing_indices(ndx,1) - 1 + n;
    peak_time = ecg_ts(virtual_peak_ndx,1);
    virtual_peak_ts(ndx,:) = [peak_time,ecg_ts(virtual_peak_ndx,2)];
    
    % extract the fitting data for peak
    window_frame = max(1,virtual_peak_ndx-peak_fitting_window):min(size(ecg_ts,1),virtual_peak_ndx+peak_fitting_window);
    virtual_peak_time = ecg_ts(window_frame,1) - ...
        ecg_ts(virtual_peak_ndx,1);
    virtual_peak_value = ecg_ts(window_frame,2);
    virtual_poly = polyfit(virtual_peak_time,virtual_peak_value,peak_fitting_order);
    virtual_polyder = polyder(virtual_poly);
    virtual_roots = roots(virtual_polyder);
    for ndx_roots = 1:length(virtual_roots)
        if abs(virtual_roots(ndx_roots)) < virtual_peak_time(end)
            true_peak_ts(ndx,1) = ecg_ts(virtual_peak_ndx,1) + virtual_roots(ndx_roots);
            true_peak_ts(ndx,2) = polyval(virtual_poly,virtual_roots(ndx_roots));
        end
    end
end
if plot_on==1
    figure(fig2),
    plot(true_peak_ts(:,1),true_peak_ts(:,2),'k*'); xlim([true_peak_ts(1,1),true_peak_ts(end,1)]);
    hold off;
end
