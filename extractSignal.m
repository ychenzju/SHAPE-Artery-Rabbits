function [signalFiltered,indices,mu,sigma] = extractSignal(timesq,signal,timeFrame,principleNum)
%% extract a segment of signal that has the most robust behavior

if nargin<4
    principleNum = 1;
end

L=round(timeFrame / 2 / ((timesq(end)-timesq(1))/(length(timesq)-1)));
Np=principleNum;
[signalDecomposed] = ssa(signal,L,Np);
signalFiltered = sum(signalDecomposed(:,1:principleNum),2);
signalVariability = nan(size(signalFiltered));
for i=1:length(timesq)
    indices = (timesq>(timesq(i)-timeFrame)) & ...
        (timesq<(timesq(i)+timeFrame));
    if sum(indices)>0 && (timesq(i)>timeFrame) && ...
            (timesq(i)<timesq(end)-timeFrame)
        signalVariability(i) = max(signalFiltered(indices)) ...
            - min(signalFiltered(indices));
    end
end
[~,indices] = min(signalVariability,[],'omitnan');
idx = indices(1);
indices = (timesq>(timesq(idx)-timeFrame)) & ...
    (timesq<(timesq(idx)+timeFrame));

mu = mean(signal(indices)); 
sigma = std(signal(indices));
