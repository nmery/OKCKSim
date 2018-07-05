function dif=gocovardm(x0libre,x,g,w,model,c,id1,id2)

%------------------------------------------------------------------------------------------------
% Function to calculate the variogram values and the differences with the experimental variogram
%------------------------------------------------------------------------------------------------

modelt=model;
ct=c;
n1=sum(sum(id1));n2=sum(sum(id2));
modelt(id1)=x0libre(1:n1);
ct(id2)=x0libre(n1+1:end);
k=covardm([0,0,0],x,modelt,ct);

gt=sum(ct)-k';
bign=99999;
dif=sum(abs(gt-g).*w)+bign*abs(sum((x0libre<0).*x0libre));