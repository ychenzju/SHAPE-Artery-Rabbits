function [A,Var_R,Residuals]=ls_AR(DAT,p)

% Dat should be a column vector

if (nargin==1),p=1;end

[N,M]=size(DAT);
if (M<1)|(M>1), 
   disp('ls_AR--------The DAT format is incorrect!');
   return;
end
  
  ls_CA=zeros(p,p);
  ls_CB=zeros(p,1);
   for ls_row=0:p-1
      for ls_col=0:p-1
         ls_a_sum=0;
         for ls_t=p:N-1
            ls_a_sum=ls_a_sum+...
               DAT(ls_t-ls_row)*DAT(ls_t-ls_col);
         end
         ls_CA(ls_row+1,ls_col+1)=ls_a_sum;
      end
      
      ls_b_sum=0;
      for ls_t=p:N-1
         ls_b_sum=ls_b_sum+...
            DAT(ls_t+1)*DAT(ls_t-ls_row);
      end
      ls_CB(ls_row+1,1)=ls_b_sum;      
   end
   ls_CA=ls_CA/N;
   ls_CB=ls_CB/N;

	ls_a=ls_CA\ls_CB;
   
   ls_x_pre=[];
   for i=p+1:N
      tsum=0;
      for k=1:p
      	tsum=tsum+ls_a(k)*DAT(i-k);
      end
      ls_x_pre=[ls_x_pre;tsum];
   end
   ls_x=DAT(p+1:N)-ls_x_pre;

   A=ls_a;
   if nargout>=2,  
      Var_R=var(ls_x);
      if nargout>=3,
         Residuals=ls_x;
      end
   end
   