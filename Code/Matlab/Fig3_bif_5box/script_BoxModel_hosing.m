%% Richard Wood box model - with hosing
clear;
close all;

%% Initial conditions

S0 = 0.035;

SN = (0.034912-S0)*100;
ST = (0.035435-S0)*100;
SS = (0.034427-S0)*100;
SIP = (0.034668-S0)*100;

H = 0;
%H = 0.384;

%% Set up for initial value problem solver

x0 = [SN;ST;SS;SIP];
tspan = [0,1000];
h = 0.1;

%% Solve the ODE

[X,t,xeq] = MyIVP(@(t,x)BoxModel_IVP_hosing(t,x,H,'FamousB_1xCO2'),x0,tspan,h);

lambda = 2.79e7; %m^6kg^-1s^-1
alpha = 0.12; %kg m^-3 C^-1
beta = 790.0;%kg m^-3 
Ts = 7.919; 
TN = 6.679;

Q = lambda.*(alpha.*(Ts-TN)+beta.*(X(1,:)./100-X(3,:)./100));

figure; plot(t,Q,'Linewidth',3);
set(gca,'FontSize',16)
%legend('SN','ST','SS','SIP');
title('Model solution: FN = 0.384');
xlabel('time (years)');
ylabel('Q');

figure; plot(X(1,:),X(2,:),'Linewidth',3);
set(gca,'FontSize',16)
%legend('SN','ST','SS','SIP');
title('Phase plane: FN = 0.384');
xlabel('SN');
ylabel('ST');


%% Bifurcation analysis- nontrivial eq 1
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.001,'h_max',0.005,'h_min',0.00001,'ItMX',2000);
prob = coco_set(prob,'corr','MX',10);
bd = coco(prob,'eq1','ode','isol','ep',... 
    @(x,p)BoxModel_coco_hosing(x,p,'FamousB_1xCO2'),xeq,{'H' 'S0'},[H S0],...  
    'H',[-0.6 0.6]);

% bd = coco_bd_read('eq1');

figure(3); hold on;
[SN1,HB1,BP1] = bifdiag1D(bd,1,'H');

%% Continuing periodic orbits
prob = coco_prob();
prob = coco_set(prob,'cont','NAdapt',1,'PtMX',50,'ItMX',[4000 0],'h0',0.01,'h_max',0.5,'h_min',0.01);
bd_po = coco(prob,'po1','ode','HB','po',... %branch, toolbox family, initial point, point type/branch type
    'eq1', HB1,... primary branch name, solution label
    {'H'},{[-0.6 0.6]}); %continuation parameter, continuation domain

figure(3); hold on;
FP1 = plotPO1D(bd_po,1,1,'po1');

% bd_po = coco_bd_read('po1');

%% Zoom for plot
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.0001,'h_max',0.0001,'h_min',0.00001,'ItMX',2000);
prob = coco_set(prob,'corr','MX',10);
bd2 = coco(prob,'eq1_zoom_2xCO2','ode','ep','ep',... 
    'eq1',19,...  
    'H',[0.16 0.4]);

% bd2 = coco_bd_read('eq1_zoom');

figure; hold on;
[SN2,HB2,BP2] = bifdiag1D(bd2,1,'H');
plotPO1D(bd_po,1,1,'po1');
box on;
set(gca,'FontSize',16)
xlim([0.2 0.23])
%% Lower equilibrium (q < 0)
% 
% [X_l,t_l,xeq_l] = MyIVP(@(t,x)BoxModel_IVP(t,x,1),x0,tspan,h);
% 
% Q_l = lambda.*(alpha.*(Ts-TN)+beta.*(X_l(1,:)-X_l(3,:)));
% 
% figure; plot(t_l,Q_l,'Linewidth',3);
% set(gca,'FontSize',16)
% %legend('SN','ST','SS','SIP');
% title('Model solution: FN = 1');
% xlabel('time (years)');
% ylabel('Q');
% 
% %% Bifurcation analysis- nontrivial eq 2
% prob = coco_prob();
% prob = coco_set(prob,'cont','h0',0.001,'h_max',0.005,'h_min',0.001,'ItMX',400);
% prob = coco_set(prob,'corr','MX',10);
% bd_l = coco(prob,'eq2','ode','isol','ep',... 
%     @BoxModel_coco_hosing,xeq_l,{'FN' 'S0'},[1 S0],...  
%     'FN',[-0.4 1]);
% 
% figure(3); hold on;
% [SN2,HB2,BP2] = bifdiag1D(bd_l,1,'FN');
% xlabel('FN')
% ylabel('SN')
% set(gca,'FontSize',20)
% xlim([-0.4 1])


%% Continuing Fold point
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.005,'h_max',0.02,'h_min',0.001,'ItMX',400);
bd_sn1 = coco(prob,'sn1','ode','SN','SN',... 
    'eq1',SN1(1),... %'eq2',HB2,...%
    {'H' 'S0'},{[-0.4 0.4] [0.005 0.07]});

prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.005,'h_max',0.02,'h_min',0.001,'ItMX',400);
bd_sn2 = coco(prob,'sn2','ode','SN','SN',... 
    'eq1',SN1(2),... %'eq2',HB2,...%
    {'H' 'S0'},{[-0.4 0.4] [0.005 0.07]});

%% Continuing Hopf point
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.005,'h_max',0.02,'h_min',0.001,'ItMX',400);
bd_hb = coco(prob,'hb1','ode','HB','HB',... 
    'eq1',HB1,... %'eq2',HB2,...%
    {'H' 'S0'},{[-0.4 0.4] [0.005 0.07]});


%% Plotting

FN_1 = coco_bd_col(bd,'H');
X_1 = coco_bd_col(bd,'x');

Q_1 = lambda.*(alpha.*(Ts-TN)+beta.*(X_1(1,:)-X_1(3,:)));

figure(5); clf; hold on;
plot(FN_1,Q_1,'b','Linewidth',3)
set(gca,'FontSize',16)
title('Bifurcation diagram: H vs Q');
xlabel('H (Sv)');
ylabel('Q (Sv)');

figure; hold on;
bifdiag1D(bd,1,'H');
xlabel('H')
ylabel('S_{N} (scaled)')
xlim([-0.5 0.5])

figure;
subplot(2,2,1); hold on;
bifdiag1D(bd,1,'H');
plotPO1D(bd_po,1,1,'po1');
xlabel('H')
ylabel('S_{N} (scaled)')
xlim([-0.4 0.4])
subplot(2,2,2); hold on;
bifdiag1D(bd,2,'H');
plotPO1D(bd_po,2,1,'po1');
xlabel('H')
ylabel('S_{T} (scaled)')
xlim([-0.4 0.4])
subplot(2,2,3); hold on;
bifdiag1D(bd,3,'H');
plotPO1D(bd_po,3,1,'po1');
xlabel('H')
ylabel('S_{S} (scaled)')
xlim([-0.4 0.4])
subplot(2,2,4); hold on;
bifdiag1D(bd,4,'H');
plotPO1D(bd_po,4,1,'po1');
xlabel('H')
ylabel('S_{IP} (scaled)')
xlim([-0.4 0.4])

%% 2 Param bifurcation diagram

br = [160 82 45]./255;
dg = [77 149 66]./225;

FN_hb = coco_bd_col(bd_hb,'FN');
S0_hb = coco_bd_col(bd_hb,'S0');

FN_sn1 = coco_bd_col(bd_sn1,'FN');
S0_sn1 = coco_bd_col(bd_sn1,'S0');

FN_sn2 = coco_bd_col(bd_sn2,'FN');
S0_sn2 = coco_bd_col(bd_sn2,'S0');

figure;
plot(FN_hb,S0_hb,'color',br);
hold on
plot(FN_sn1,S0_sn1,'color',dg);
plot(FN_sn2,S0_sn2,'color',dg);
xlabel('FN')
ylabel('S0')
title('Brown- Hopf, Green - Fold')