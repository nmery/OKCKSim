function [be,bt,osre,osrt,neg,ke,be_cr,bt_cr,osre_cr,osrt_cr,ke_cr,stat]=test_ck(n,modelg,cg,model,c,l2,m,tab)

%-------------------------------------------------------------------------------------------------------
% Function to perform the simulations, compute KNA statistics, recovery functions and other stats.
%-------------------------------------------------------------------------------------------------------
%
% USE
%      grille       : function to create a block discretization
%      covardm      : function to calculate the covariance (Marcotte, 1991)
%
% INPUT
%      n            : number of simulated block
%      modelg       : initial proposed covariance model.
%      cg           : initial coefficient matrix of the coregionalization model.
%      model        : covariance model.
%      c            : coefficient matrix of the coregionalization model.
%      l2           : sampling density.
%      m            : vector with the number of neighbors.
%      tab          : database.
%
% OUTPUT
%      be           : experimental slope of regression (OK).
%      bt           : theoretical slope of regression (OK).
%      osre         : experimental over smoothing ratio (OK). 
%      osrt         : theoretical over smoothing ratio (OK).
%      neg          : proportion of negative weights (OK).
%      ke           : kriging efficiency (OK).
%      be_cr        : experimental slope of regression (CK).
%      bt_cr        : theoretical slope of regression (CK).
%      osre_cr      : experimental over smoothing ratio (CK).
%      osrt_cr      : theoretical over smoothing ratio (CK).
%      ke_cr        : kriging efficiency (CK).
%      stat         : main statistics.
%-------------------------------------------------------------------------------------------------------

% Definition of constants and parameters       
nn=125+m;
zv=zeros(n,1);
ze=zeros(n,1);
ze_cr=ze;
bt1=zeros(n,1);
bt1_cr=bt1;
osrt=zeros(n,1);
osrt_cr=osrt;
s=zeros(n,1);
neg=zeros(n,1);
KE=zeros(n,1);
KE_cr=KE;
index=zeros(n,1);
tab=[-7 0;tab;7 1.5*max(tab(:,2))];
nbimag=0;

% Block discretization
x0=grille(-2,2,1,-2,2,1,-2,2,1);

% Compute block covariance
k0=covardm(x0,x0,model,c); 
sv=mean(k0(:));   

% Simulation of points and compute kriging (ordinary and constrained)
for i=1:n
    x=rand(m,3)*l2-l2/2;                              % localize randomly points around the block
    kx=covardm(x,x,model,c);                          % point covariance
    kxx0=covardm(x,x0,model,c);                       % point-block covariance
    k=[kx kxx0;kxx0' k0]; 
    kg=[kx ones(m,1);ones(1,m) 0];                    % left hand kriging matrix
    kd=[mean(kxx0,2);1];                              % right hand kriging vector
    l=kg\kd;                                          % kriging weights
    ll=l(1:m);                                        % kriging weights without Lagrange multiplier
    varze=ll'*kx*ll;                                  % theoretical variance of the estimate Var(Z*)
    t=diag(l'*kd);
    s(i)=sv-t';
    
        % Ordinary kriging
         bt1(i)=100*ll'*kd(1:m)/varze;                     % theoretical slope of regression
         osrt(i)=100*((sv-varze)/sv);                      % theoretical oversmoothing ratio
         neg_n=find(ll<0);                                 % negative weights
         negg_n=ll(neg_n);
         count_neg=length(negg_n);
         neg(i)=100*count_neg./length(ll);
         KE(i)=100.*(sv-s(i))./sv;                         % kriging efficiency
    
        % Constrained kriging
         ki=inv(kx);                                       % inverse of variance-covariance matrix
         ski=sum(ki(:));                                   % sum of all elements in ki
         l=ki*kd(1:m);
         k0l=kd(1:m)'*l;
         mu2=sqrt(l'*kd(1:m)*ski-sum(l)^2)/sqrt(ski*sv-1); % Eq. 3.13 (Cressie, 1993b)
         mu1=-mu2/ski+sum(l)/ski;                          % Eq. 3.12 (Cressie, 1993b)
         lc=1/mu2*kd(1:m)'*ki-mu1/mu2*sum(ki);             % Eq. 3.11 (Cressie, 1993b)
         lc=lc';
         if sum(abs(imag(lc)))>0                           % Case for complex number in CK
            s_cr(i)=NaN;
            bt1_cr(i)=NaN;
            osrt_cr(i)=NaN;
            KE_cr(i)=NaN;
            index(i)=0;
            varze_cr=NaN;   
            nbimag=nbimag+1;
            lc=ll;                                         % when complex number appears using CK take ordinary kriging estimate   
         else
            varze_cr=lc'*kx*lc;                               % theoretical variance of the estimate Var(Z*)
            s_cr(i)=sv+varze_cr-2*lc'*kd(1:m);
            bt1_cr(i)=100*lc'*kd(1:m)/varze_cr;               % theoretical slope of regression
            osrt_cr(i)=100*((sv-varze_cr)/sv);                % theoretical over-smoothing ratio
            KE_cr(i)=100.*(sv-s_cr(i))./sv;                   % kriging efficiency
            index(i)=i;
         end
         
    % Simulation of block and sample values
    kx=covardm(x,x,modelg,cg);
    kxx0=covardm(x,x0,modelg,cg);
    k0=covardm(x0,x0,modelg,cg);
    k=[kx kxx0;kxx0' k0];
    u=chol(k);                                             
    y=u'*randn(nn,1);                                  % block discretization points and sample are simulated by Cholesky
    z=interp1(tab(:,1),tab(:,2),y,'linear','extrap');
    zv(i)=mean(z(m+1:end));                            % true block value
    ze_cr(i)=lc'*z(1:m);                               % constrained kriging estimate
    ze(i)=ll'*z(1:m);                                  % ordinary kriging estimate 
end

 disp(['Computing KNA criteria with ',num2str(m),' neighbors ']);
 
% Remove negative values
ze_cr(ze_cr<0)=0;
ze(ze<0)=0;
index(index==0)=[];

% Mean of estimate values [true OK CK]
[mean(zv) mean(ze) mean(ze_cr)];

% Compute experimental and mean of theoretical KNA criteria
be=regres(zv,ze,1);                                    % experimental slope of regression (OK)
be=be(1)*100;
bt=nanmean(bt1);                                       % mean of theoretical slope of regression (OK)
osrt=nanmean(osrt);                                    % mean of theoretical over-smoothing ratio (OK)
neg=nanmean(neg);                                      % mean of proportion of negative values (OK)
ke=nanmean(KE);                                        % mean of kriging efficiency (OK)
osre=100*((var(zv)-var(ze))/var(zv));                  % experimental over-smoothing ratio (OK)
be_cr=regres(zv(index),ze_cr(index),1);                % experimental slope of regression (CK)
be_cr=be_cr(1)*100;
if be_cr==0
    be_cr=NaN;
else
    be_cr=be_cr;
end
bt_cr=nanmean(bt1_cr);                                 % mean of theoretical slope of regression (CK)
osrt_cr=nanmean(osrt_cr);                              % mean of theoretical over-smoothing ratio (CK)  
ke_cr=nanmean(KE_cr);                                  % mean of kriging efficiency (CK)
osre_cr=100*((var(zv(index))-var(ze_cr(index)))/var(zv(index))); % experimental over-smoothing ratio (CK)

% Compute recovery functions (tonnage, ore grade, conventional profit)
p=[0.025:0.05:1]';
c=quantile(zv,p);
for i=1:length(c)
    tc(i)=mean(zv>c(i));
    tce(i)=mean(ze>c(i));
    tc_cr(i)=mean(ze_cr>c(i));
    
    mc(i)=mean(zv(zv>c(i)));
    mce(i)=mean(ze(ze>c(i)));
    mc_cr(i)=mean(ze_cr(ze_cr>c(i)));
    
    pc(i)=tc(i)*(mc(i)-c(i));
    pce(i)=tce(i)*(mce(i)-c(i));
    pc_cr(i)=tc_cr(i)*(mc_cr(i)-c(i));
end

% Statistics
stat.propneg=nanmean(ze<0);
stat.propneg_cr=nanmean(ze_cr<0);
stat.prop_imag=nbimag/n;
stat.n=n;
stat.nbneg=sum(ze<0);
stat.nbneg_cr=sum(ze_cr<0);
stat.nbimag=nbimag;
stat.rmse=sqrt(nanmean((zv-ze).^2));
stat.rmse_cr=sqrt(nanmean((zv-ze_cr).^2));
t=corrcoef(zv,ze);
stat.corr=t(1,2);
t=corrcoef(zv,ze_cr);
stat.corr_cr=t(1,2);
stat.tc=tc;
stat.mc=mc;
stat.pc=pc;
stat.tce=tce;
stat.mce=mce;
stat.pce=pce;
stat.tc_cr=tc_cr;
stat.mc_cr=mc_cr;
stat.pc_cr=pc_cr;

% Plot recovery functions

disp(['Plot recovery function considering ',num2str(m),' neighbors ']);
% Tonnage
figure(100+m);clf
plot(c,tc,'-r',c,tce,'-b',c,tc_cr,'-k','linewidth',2)
title(['Neighbors ',num2str(m)],'fontsize',15)
grid on
h=legend('Real','OK','CK');
set(h,'Interpreter','Latex','Location','northeast','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
box off
xlabel('Cut-off','Interpreter','Latex','FontSize', 15)
ylabel('Tonnage','Interpreter','Latex','FontSize', 15)
title(['Curve tonnage - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 15)

% Ore grade
figure(101+m)
plot(c,mc,'-r',c,mce,'-b',c,mc_cr,'-k','linewidth',2)
h=legend('Real','OK','CK');
set(h,'Interpreter','Latex','Location','southeast','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title(['Curve ore grade - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 15)
grid on
box off
xlabel('Cut-off','Interpreter','Latex','FontSize', 15)
ylabel('Ore grade','Interpreter','Latex','FontSize', 15)

% Conventional profit
figure(102+m)
plot(c,pc,'-r',c,pce,'-b',c,pc_cr,'-k','linewidth',2)
h=legend('Real','OK','CK');
set(h,'Interpreter','Latex','Location','northeast','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title(['Curve Conv. Profit - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 15)
grid on
box off
xlabel('Cut-off','Interpreter','Latex','FontSize', 15)
ylabel('Conventional profit','Interpreter','Latex','FontSize', 15)

drawnow