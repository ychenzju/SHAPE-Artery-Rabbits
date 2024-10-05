function [signalFiltered]=ssa(signal,windowLen,Np)
%========================================================================
% signal 原始信号
% windowLen 窗口长度
% signalFiltered   重构时间序列
%=========================================================================

if nargin<2
    windowLen = round(N/2) - 1;
elseif nargin<3
    Np = round(windowLen / 2);
end

% Step1 : 建立轨迹矩阵
N=length(signal);
if windowLen>N/2
    windowLen=N-windowLen;
end
K=N-windowLen+1;
X=zeros(windowLen,K);
for i=1:K
    X(1:windowLen,i)=signal(i:windowLen+i-1);
end

% Step 2/3: PCA
[P,T,R] = pca_pls(X',Np);
%P=U(:,I)
%T=V(:,I)

signalFiltered = zeros(N,Np);
for i=1:Np
    rca = P(:,i) * T(:,i)';
    
    % Step 4: 对交平均化重构信号
    y=zeros(N,1);
    Lp=min(windowLen,K);
    Kp=max(windowLen,K);
    %重构 1~Lp-1
    for k=0:Lp-2
        for m=1:k+1
            y(k+1)=y(k+1)+(1/(k+1))*rca(m,k-m+2);
        end
    end
    
    %重构 Lp~Kp
    for k=Lp-1:Kp-1
        for m=1:Lp
            y(k+1)=y(k+1)+(1/(Lp))*rca(m,k-m+2);
        end
    end
    
    %重构 Kp+1~N
    for k=Kp:N-1
        for m=k-Kp+2:N-Kp+1
            y(k+1)=y(k+1)+(1/(N-k))*rca(m,k-m+2);
        end
    end
    
    signalFiltered(:,i) = y;
end

end