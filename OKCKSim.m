%-------------------------------------------------------------------------------------------------------------------------------
% Main routine to drive the computations of KNA and recovery functions using ordinary kriging (OK) and constrained kriging (CK)
%-------------------------------------------------------------------------------------------------------------------------------
%
% USE
%      do_anamor    : routine to create a table associating Gaussian values to real values.
%      transcov     : function to simulate and fit the real variogram from the Gaussian variogram.
%      test_ck      : function to perform the simulations, compute KNA statistics, recovery functions and other stats.
%
% USER DEFINITIONS
%      SD_vol       : Volume per unit sample 100, 500 and 2500 for respectively HSD, Base Case and HSD (line 18).
%      cas          : To perform the analysis with the three datasets 'cas' (line 19) have to be changed. 1 for lognormal
%                     distribution, 2 for reverse lognormal distribution and 3 for bimodal distribution.
%-------------------------------------------------------------------------------------------------------------------------------

% Definition of parameters
seed=9153;              % random number initialisation
SD_vol=500;             % choice of sampling density 100, 500 and 2500 for HSD, BC and LSD in the paper, average volume per unit sample.
cas=1;                  % type of distribution 1: lognormal, 2 reverse lognormal, 3- bimodal
vn=[5 8 10 20 30 50];   % vector with the number of neighbors
n=15000;                % number of simulated blocks for each case (each block is simulated independently)

rng('default')
rng(seed);

% Create datasets
disp('Create datasets')
do_anamor;

% Write the number of the case for each distribution (1: lognormal, 2: reverse lognormal and 3: bimodal)
disp('Fitting variogram of transformed variables') 
switch cas
    case 1                                      % Lognormal distribution
        modelg=[1 1;4 20];cg=[0.3;0.7];
        yz=tab(:,:,1);
        [model,c]=transcov(modelg,cg,yz);
        
    case 2                                      % Reverse lognormal distribution
        modelg=[1 1;4 20];cg=[0.3;0.7];
        yz=tab(:,:,2);
        [model,c]=transcov(modelg,cg,yz);
        
    case 3                                      % Bimodal distribution
        modelg=[1 1;4 20];cg=[0.3;0.7];
        yz=tab(:,:,3);
        [model,c]=transcov(modelg,cg,yz);
end

% Run considering different number of neighbors (vn)
for ii=1:length(vn)
    m=vn(ii);
    disp(' ');
    disp(['Computing OK and CK estimates with ',num2str(m),' neighbors (case ',num2str(ii),' of ', num2str(length(vn)),')']);
    l2=(SD_vol*m)^(1/3);                           % sampling density (modify 500 for HSD and LSD)
    [be,bt,osre,osrt,neg,ke,be_cr,bt_cr,osre_cr,osrt_cr,ke_cr,stat]=test_ck(n,modelg,cg,model,c,l2,m,yz);
    statall{ii}=stat;                           % save statistics of ii neighbors
    vbe(ii)=be;                                 % save experimental slope of regression (OK) of ii neighbors
    vbt(ii)=bt;                                 % save theoretical slope of regression (OK) of ii neighbors
    vosre(ii)=osre;                             % save experimental over-smoothing ratio (OK) of ii neighbors
    vosrt(ii)=osrt;                             % save theoretical over-smoothing ratio (OK) of ii neighbors
    vneg(ii)=neg;                               % save negative weights (OK) of ii neighbors
    vke(ii)=ke;                                 % save kriging efficiency (OK) of ii neighbors
    vbe_cr(ii)=be_cr;                           % save experimental slope of regression (CK) of ii neighbors
    vbt_cr(ii)=bt_cr;                           % save theoretical slope of regression (CK) of ii neighbors
    vosre_cr(ii)=osre_cr;                       % save experimental over-smoothing ratio (CK) of ii neighbors
    vosrt_cr(ii)=osrt_cr;                       % save theoretical over-smoothing ratio (CK) of ii neighbors
    vke_cr(ii)=ke_cr;                           % save kriging efficiency (CK) of ii neighbors
end

% Plot KNA criteria vs number of neighbors

% Ordinary kriging results
disp(' ');
disp('Plot KNA criteria using OK');
figure(10*cas);clf
plot(vn,vbt,'o-b','linewidth',2);
hold on
plot(vn,vbe,'x--b','linewidth',2);
plot(vn,vosrt,'o-r',vn,vosre,'x--r','linewidth',2);
plot(vn,vke,'-k','linewidth',2)
axis([0 55 -40 110])
grid on
hold off
h=legend('$SR_{th}$','$SR_{exp}$','$OSR_{th}$','$OSR_{exp}$','$KE$');
set(h,'Interpreter','Latex','Location','southeast');
set(gca,'TickLabelInterpreter','latex')
box off
xlabel('Number of neighbors','Interpreter','Latex','FontSize', 12)
ylabel('$\%$','Interpreter','Latex','FontSize', 12)

% Constrained kriging results
disp('Plot KNA criteria using CK');

figure(10*cas+1);clf;
plot(vn,vbt_cr,'o-b','linewidth',2);
hold on
plot(vn,vbe_cr,'x--b','linewidth',2);
plot(vn,vosrt_cr,'o-r',vn,vosre_cr,'x--r','linewidth',2);
plot(vn,vke_cr,'-k','linewidth',2)
axis([0 55 -20 110])
grid on
hold off
h=legend('$SR_{th}$','$SR_{exp}$','$OSR_{th}$','$OSR_{exp}$','$KE$');
set(h,'Interpreter','Latex','Location','southeast');
set(gca,'TickLabelInterpreter','latex')
box off
xlabel('Number of neighbors','Interpreter','Latex','FontSize', 12)
ylabel('$\%$','Interpreter','Latex','FontSize', 12)

disp(' ');
disp('...Done ...');