function [modela,ca,dif]=ajustermulti(gexp,model,c,modelcode,ccode,vdir,vreg,expoh,hmin);

%--------------------------------------------------------------------------------------------------
% Function to simultaneously adjust a model in multiple directions (Marcotte, 2002).
%--------------------------------------------------------------------------------------------------
%
% INPUT
%      gexp         : experimental variogram  in the following order [h,n(h),g(h)].
%      model        : variogram model.
%      c            : coefficient matrix of the coregionalization model.
%      modelcode    : 0 for free or 1 for imposed adjustment of the model.
%      ccode        : 0 for free or 1 for imposed adjustment of the c for the different structures.
%      vdir         : directions used to calculate the experimental variogram.
%      vreg         : vector containing the regularizations used to calculate the exp variogram.
%      expoh        : exponent.
%      hmin         : minimun distance.
%
%
% OUTPUT
%      modela       : fitted covariance model.
%      ca           : fitted coefficient matrix of the coregionalization model.
%      dif          : adjusted statistics.
%---------------------------------------------------------------------------------------------------

% Definition of constants
[n,p]=size(vdir);
if p==1
    vdir=[vdir zeros(n,1)];
end
u=poletocart(vdir);
ndir=size(vdir,1);
x=[];
g=[];
w=[];

%Weight by number of pairs and distance
for idir=1:ndir
    id=gexp(:,2,idir)>0;
    x=[x; gexp(id,1,idir)*u(idir,:)];
    g=[g; gexp(id,3,idir)];
    w=[w;gexp(id,2,idir).*ones(sum(id),1)./(hmin+gexp(id,1,idir)).^expoh];
end

[n,p]=size(model);
deuxd=0;

% Check if 2D or 3D anisotropy
if p==4
    id=model(:,2:3)==0;
    model(id)=0.0001;
    deuxd=1;
    model=[model(:,1:3),model(:,3),zeros(n,2),model(:,4)];
    modelcode=[modelcode(:,1:3),ones(n,3),modelcode(:,4)];
elseif p==7
    id=model(:,2:4)==0;
    model(id)=0.0001;
end

id=model(:,1)==1;
model(id,2:end)=1;
modelcode(id,2:end)=1;
modelcode(:,1)=1;

options=optimset('TolX',0.00001,'TolFun',0.00001,'MaxFunEvals',13500);

id1=modelcode==0;p1=sum(sum(id1));
x0libre=reshape(model(id1),p1,1);
id2=ccode==0;p2=sum(sum(id2));
x0libre=[x0libre;c(id2)];

if ~isempty(x0libre)
    [t1,dif,flag]=fminsearch('gocovardm',x0libre,options,x,g,w,model,c,id1,id2);
    if flag==0
        display('Stop on number of iterations')
    else
        display('Stop at a local minimum to simultaneously adjust a covariance model')
    end
    if p1>0
        t11=t1(1:p1);
        model(id1)=t11;
    end
    if p2>0
        t21=t1(p1+1:p1+p2);
        c(id2)=t21;
    end
else
    dif=gocovardm(x0libre,x,g,w,model,c,id1,id2);
    modela=model;
    ca=c;
end

if deuxd==1
    modela=model(:,[1:3 7]);
else
    modela=model;
end
ca=c;

%Plot the variograms
ndir=size(vdir,1);
nplot=ceil(sqrt(ndir));
if nplot*nplot==ndir
    nplot=nplot+1;
end
dh=[];hmax=[];gmax=[];
for idir=1:ndir
    id=gexp(:,2,idir)>0;
    dh=[dh;min(gexp(id,1,idir))];
    hmax=[hmax;max(gexp(id,1,idir))];
    gmax=[gmax;max(gexp(id,3,idir))];
end
dh=0.1*min(dh);
hmax=1.2*max(hmax);
gmax=1.2*max(gmax);
h=[dh:dh:hmax]';
u=poletocart(vdir);

for idir=1:ndir
    figure(4);clf;
    nn=find(gexp(:,2,idir)>0);
    id=gexp(:,2,idir)>0;
    x=h*u(idir,:);
    k=covardm(x,[0,0,0],model,c);
    gt=sum(c)-k;
    set(gca,'TickLabelInterpreter','latex')
    box off
    plot(gexp(id,1,idir),gexp(id,3,idir),'+','markersize',5);hold on
    plot(h,gt);hold off
    axis([0 hmax, 0, gmax]);
    set(gca,'TickLabelInterpreter','latex','FontSize',15)
    box off
    title('Variogram dataset 1','Interpreter','Latex','FontSize',15)
    xlabel('Lag distance','Interpreter','Latex','FontSize', 15)
    ylabel('$\gamma (h)$','Interpreter','Latex','FontSize', 15)
    h=legend('Experimetal variogram','Fitted variogram');
    set(h,'Interpreter','Latex','Location','southeast','FontSize',15);
end

if isempty(findobj(gcf,'label','Affichage'))
    t=uimenu(gcf,'label','Affichage');
    t1=uimenu(t,'label','nombre de paires','checked','on','tag','apaire',...
        'callback','apaire');
    t2=uimenu(t,'label','modèle','checked','on','tag','amodele',...
        'callback','amodele2');
end

function x=poletocart(pole)
deg2rad=pi/180;
pole=pole*deg2rad;

x(:,1)=cos(pole(:,2)).*sin(pole(:,1));
x(:,2)=cos(pole(:,2)).*cos(pole(:,1));
x(:,3)=-sin(pole(:,2));