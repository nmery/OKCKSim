function [model,c]=transcov(modelg,cg,yz)

%-------------------------------------------------------------------------------------------------------
% Function to simulate and fit the real variogram from the Gaussian variogram (Marcotte, 2002).
%-------------------------------------------------------------------------------------------------------
%
% USE 
%      covardm      : function to calculate the covariance (Marcotte, 1991)
%      ajustermulti : function to simultaneously fit a model in multiple directions (Marcotte, 2002)
%
% INPUT
%      modelg       : initial proposed model structure. The first column a code for the model type 
%                     (1: nugget effect; 2: exponential; 3: Gaussian; 4: spherical), the following
%                     columns give the ranges along the different coordinates.
%      cg           : coefficient matrix of the coregionalization model.
%      yz           : values of each dataset (read from tab.mat).
%
% OUTPUT
%      model        : fitted covariance model. 
%      c            : fitted coefficient matrix of the coregionalization model.
%-------------------------------------------------------------------------------------------------------

%Parameters
seed=915;
rng(seed);
yz=[-7 0;yz;7 1.5*max(yz(:,2))];
hmax=1.2*max(modelg(:,2));
dh=hmax/100;
h=[0, 0.00001:dh:hmax]';
n=100000;

%Compute the covariance
k=covardm(0,h,modelg,cg);
y=randn(n,1);
z=interp1(yz(:,1),yz(:,2),y,'linear','extrap');
kt=zeros(size(k));
kt(1)=var(z);

for i=2:length(k)
    kk=[k(1) k(i);k(i) k(1)];
    u=chol(kk);
    yy=u'*randn(2,n);
    t1=interp1(yz(:,1),yz(:,2),yy(1,:)');
    t2=interp1(yz(:,1),yz(:,2),yy(2,:)');
    tt=cov(t1,t2);
    kt(i)=tt(1,2);
end

tt=kt(1)-kt(2:end);
gexp=[h(2:end) ones(length(h)-1,1) tt(:)];
[model,c,dif]=ajustermulti(gexp,[1 1;2 hmax/6;4 hmax/2],kt(1)*[1/3;1/3;1/3],[1 1;1 0;1 0],[0;0;0],[0 0],1,1,1);