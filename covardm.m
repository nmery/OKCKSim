function [k]=covardm(x,x0,model,c)

%--------------------------------------------------------------------------------------------------------------------------
% Function to calculate the covariance with specified models and with cokri (Marcotte, 1991)
%--------------------------------------------------------------------------------------------------------------------------
%
% USE 
%      trans        : generate rotated and reduced coordinates following specifications described in model (Marcotte, 1997)
%      
% INPUT
%      x            : data matrix
%      x0           : matrix of coordinates 
%      model        : proposed model structure. The first column a code for the model type 
%                     (1: nugget effect; 2: exponential; 3: Gaussian; 4: spherical), the following
%                     columns give the ranges along the different coordinates.
%      c            : coefficient matrix of the coregionalization model.
%
% OUTPUT
%      k            : covariance matrix  
%-------------------------------------------------------------------------------------------------------------------------

% Covariance models
Gam=['h==0                                                               '; %nugget
     'exp(-h)                                                            '; %exponential
     'exp(-(h).^2)                                                       '; %gaussian
     '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)                              '; %spherical
     '1-h                                                                '; %linear
     '1-(7*min(h,1).^2-8.75*min(h,1).^3+3.5*min(h,1).^5-0.75*min(h,1).^7)'; %modele cubique
     '(h.^2).*log(max(h,eps))                                            '; %spline plaque mince
     '(h.^2+1).^(-0.5)                                                   '; %modèle gravimétrique (Cauchy avec b=0.5)
     '(h.^2+1).^(-1.5)                                                   '; %modele magnétique (Cauchy avec b=1.5) 
     'sin(max(eps,h*2*pi))./max(eps,h*2*pi)                              '; %effet de trou sinusoidal
     'cos(h*2*pi)                                                        '; %effet de trou cosinusoidal
     '1-h.^2./(1+h.^2)                                                   '];%christakos 1984

% Definition of constants
k=[];
[n1,d]=size(x);
[n2,d]=size(x0);
[rp,p]=size(c);
r=rp/p;
cx=[x(:,1:d);x0];
nm=size(model,2);

% Do not allow scopes of 0 in input to avoid divisions by 0
if nm>2
    model(:,2:1+d)=max(model(:,2:1+d),100*eps);
else
    model(:,2)=max(model(:,2),100*eps);
end

% Compute covariance
k=zeros(n1*p,n2*p);

for i=1:r
    [t1]=trans(x(:,1:d),model,i);
    [t2]=trans(x0,model,i);
    h=0;
    for id=1:d
        h=h+(t1(:,id)*ones(1,n2)-ones(n1,1)*t2(:,id)').^2;
    end
    h=sqrt(h);
    ji=(i-1)*p+1; js=i*p ;
    g=eval(Gam(model(i,1),:));
    k=k+kron(g,c(ji:js,:));
end