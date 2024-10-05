function [Beta,F,dim]=pls_modeling0(X,Y,mode,var_limit)

% Y and X must be with mean 0.

% Y=X*Beta+F
% mode 0: Utilizing CROSS_VALIDITY method to determine the model dimension
% mode 1: Utilizing VARIANCE_CONTRIBUTION_RATE method to determine the 
%         model dimension, and the rate limit is set as var_limit
%         which is preset as 0.85 defaultly.
% pls_modeling(X,Y)====pls_modeling(X,Y,0)====pls_modeling(X,Y,k),k!=1
% pls_modeling(X,Y,1,0.85)

if (nargin==2),
   mode=0;
elseif (nargin==3),
   if (mode==1),
      var_limit=0.85;
   else
      mode=0;%other values are viewd as mode 0 
   end
elseif (nargin==4),
   if (mode==0)
      disp('The mode must be 1 if nargin>=4');
      mode=1;
   else
      mode=1;%other values are viewd as mode 1 
      if (var_limit>1)|(var_limit<0) 
         disp('The precision of regression is expressed by ');
         disp('the rate between variance fitted and oringinal variance');
         disp('As a result, the var_limit must be a value between 0 and 1');
         disp('In this situation, var_limit is set to 0.85 defaultly')
      end
   end
end

[n1,m1]=size(X);
[n2,m2]=size(Y);
if (n1<n2)|(n1>n2),
   disp('The length of X must be equal to the one of Y!');
   return;
end
n=n1;

if (mode==0)
   
   part=5;
   partition(1)=1;
   for i=2:part
      partition(i)=round((i-1)*n/part);
   end
   partition(part+1)=n;
   
   press_total=10000000;
   press_total_last=inf;
   
   h=0;
   while ((press_total<press_total_last)&(h<m1))
      h=h+1;
      sumOfPress=0;
      for j=1:part
         if (j==1),
            XX=X(partition(2):partition(part+1),:);
            YY=Y(partition(2):partition(part+1),:);
         elseif (j<part)
            XX=[X(partition(1):partition(j),:);...
                  X(partition(j+1):partition(part+1),:)];
            YY=[Y(partition(1):partition(j),:);...
                  Y(partition(j+1):partition(part+1),:)];
         else
            XX=X(1:partition(part),:);
            YY=Y(1:partition(part),:);
         end
         
         [W,P,Q,B,RYY,RXX]=pls(XX,YY,h);
         XXX=X(partition(j):partition(j+1),:);
         YYY=Y(partition(j):partition(j+1),:);
         %Beta=W*inv(P'*W)*Q';
         %Beta=P*B*Q';
         %predict_YYY=XXX*Beta;
         predict_YYY=pls_predictor(XXX,P,B,Q);
         
         sumOfPress=sumOfPress+press(YYY-predict_YYY);
      end
      
      press_total_last=press_total;
      press_total=sumOfPress;
      press_with_h=[h,press_total]
   end
   
   if (h<2)|(press_total<press_total_last)
      [W,P,Q,B,RY,RX]=pls(X,Y,h);
      if (nargout==3), dim=h;end
   else
      [W,P,Q,B,RY,RX]=pls(X,Y,h-1);
      if (nargout==3), dim=h-1;end
   end
   
   Beta=W*inv(P'*W)*Q';   
   %Beta=P*B*Q';
   
   F=RY;
   PLS_T_H=var(X*P)
   %sum(PLS_T_H(1:3))/sum(PLS_T_H)
   
else
   
   [W,P,Q,B,RY,RX]=pls(X,Y,m1);
   PLS_T_H=var(X*P)
   
   dimension=0;
   var_limit
   for dd=1:m1
      contribution_rate=sum(PLS_T_H(1:dd))/sum(PLS_T_H);
      if (contribution_rate>var_limit),
         dimension=dd;
	      break;
      end
   end
   
   if (dimension>0)
      [W,P,Q,B,RY,RX]=pls(X,Y,m1);
   end
   %dimension
   %contribution_rate=sum(PLS_T_H(1:dimension))/sum(PLS_T_H)

   Beta=W*inv(P'*W)*Q';
   F=RY;
   
   if (nargout==3)
      dim=dimension;
   end
   
end




      
   

      

