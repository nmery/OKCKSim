function [b,yc,e,sc,r2,s2b,rb]=regres(y,x,ord)

%--------------------------------------------------------------------------------------------------------------------------
% Function to compute the regression of x over y 
%--------------------------------------------------------------------------------------------------------------------------
%      
% INPUT
%      x            : independet variable
%      y            : dependent variable 
%      ord          : 0 for model without constant; 1 for a constant in the model
%
% OUTPUT
%      b            : vector of model coefficients   
%      yc           : vector with predicted values  
%      e            : vector with errors (yc-y) 
%      sc           : ANOVA table  
%      r2           : multiple correlation coefficient
%      s2b          : variance-covariance matrix of the coefficients b 
%      rb           : matrix of correlations between the coefficients b 
%-------------------------------------------------------------------------------------------------------------------------

% Definition of parameters       
[n,p]=size(x);
if ord==1
    x=[x,ones(n,1)];
end
b=x\y;			                         	% coefficients of regression
yc=x*b;                                     % predicted values
e=y-yc;                                     % errors
nl=(ord+1)*3;                               % number of lines 
sc=ones(nl,3)*nan;
sc(1,1)=y'*y;   sc(1,3)=n;                  % SCT and d.l.
sc(2,1)=yc'*yc; sc(2,3)=p+ord;              % SCR and d.l.
sc(3,1)=e'*e;   sc(3,3)=sc(1,3)-sc(2,3);    % SCE and d.l.
if ord==1
    sc(4,1)=n*mean(y)^2; sc(4,3)=1;         % SCM and d.l.
    sc(5,1)=sc(1,1)-sc(4,1) ; sc(5,3)=n-1;  % SCTm and d.l.
    sc(6,1)=sc(2,1)-sc(4,1) ; sc(6,3)=p;    % SCRm and d.l.
end
sc(:,2)=sc(:,1)./max(sc(:,3),1);            % CM
if ord==1
    r2=sc(6,1)/sc(5,1);                     % R2
else
    r2=1-sc(3,1)/(sc(1,1)-n*mean(y)^2);     % new way to calculate R2 without constant
end
s2b=sc(3,2)*inv(x'*x);                      % Variance-cov of b coeff. 
t=sqrt(diag(s2b)); t=inv(diag(t));
rb=t*s2b*t;                                 % correlation of b coeff. 
end