%--------------------------------------------------------------------------
% Routine to create a table associating Gaussian values to real values
%--------------------------------------------------------------------------
%
% OUTPUT
%      tab.mat      : table with 3 different anamorphosis
%--------------------------------------------------------------------------

% Create three anamorphosis
clear tab;
N=100000;
p=[0.000001:0.000001:0.00001 0.00002:0.00001:0.0001 0.0002:0.0001:0.001 0.002:0.001:0.01 0.02:0.01:0.99 0.991:0.001:0.999 0.9991:0.0001:0.9999 0.99991:0.00001:0.99999 0.999991:0.000001:0.999999]';
y=norminv(p,0,1);
tab(:,1,1)=y;
tab(:,1,2)=y;
tab(:,1,3)=y;

% Lognormal distribution
z=exp(randn(N,1));
m=mean(z);
s=std(z);
sk=skewness(z);
tab(:,2,1)=quantile(z,p);
zz=z(z<=10);
figure(1);clf
hist(zz,100);
box off;
set(gca,'TickLabelInterpreter','latex','FontSize',15);
str={['$\mu$=',num2str(sprintf('%.1f',m))];['$\sigma$=',num2str(sprintf('%.1f',s))];['sk=',num2str(sprintf('%.1f',sk))]};
text(8.1,6000,str,'Interpreter','latex','FontSize', 15,'EdgeColor','black');
xlabel('Data values','Interpreter','Latex','FontSize', 15);
ylabel('Number of data','Interpreter','Latex','FontSize', 15);

% Reversed lognormal distribution
z=betarnd(3*ones(N,1),1.2);
z=z/std(z)*s;
z=z-min(z);
m=mean(z);
s=std(z);
sk=skewness(z);
tab(:,2,2)=quantile(z,p);
figure(2); clf
hist(z,100);
xlim([0 10]);
box off;
set(gca,'TickLabelInterpreter','latex','FontSize',15);
str={['$\mu$=',num2str(sprintf('%.1f',m))];['$\sigma$=',num2str(sprintf('%.1f',s))];['sk=',num2str(sprintf('%.1f',sk))]};
text(0.32,2150,str,'Interpreter','latex','FontSize', 15,'EdgeColor','black');
xlabel('Data values','Interpreter','Latex','FontSize', 15);
ylabel('Number of data','Interpreter','Latex','FontSize', 15);

% Bimodal distribution
z=[exp(randn(N,1));1.1*randn(N/2,1)+7];
z=max(0,z);
z=min(50,z);
zz=z(z<=10);
m=mean(z);
s=std(z);
sk=skewness(z);
tab(:,2,3)=quantile(z,p);
figure(3);clf
hist(zz,100);
box off;
set(gca,'TickLabelInterpreter','latex','FontSize',15);
str={['$\mu$=',num2str(sprintf('%.1f',m))];['$\sigma$=',num2str(sprintf('%.1f',s))];['sk=',num2str(sprintf('%.1f',sk))]};
text(8.3,6000,str,'Interpreter','latex','FontSize', 15,'EdgeColor','black');
xlabel('Data values','Interpreter','Latex','FontSize', 15);
ylabel('Number of data','Interpreter','Latex','FontSize', 15);

% Save tab file
save tab.mat   