%% Richard Wood box model - 2D with hosing
clear;
close all;

%% Initial conditions

S0 = 0.035;

SN = (0.034912-S0)*100;
ST = (0.035435-S0)*100;
SS = (0.034427-S0).*100;

H = 0;
%gamma = 0.39;
gamma = 0.36; % 2 x CO2

%% Set up for initial value problem solver

x0 = [SN;ST];
tspan = [0,2000];
h = 0.1;

%% Solve the ODE

[X,t,xeq] = MyIVP(@(t,x)BoxModel_2DH_IVP(t,x,H,'FamousB_2xCO2'),x0,tspan,h);

alpha = 0.12; %kg m^-3 C^-1
beta = 790.0;%kg m^-3 
Ts = 7.919; 
T0 = 2.65; %C
lambda = 2.79e7; %m^6kg^-1s^-1
mu = 5.5e-8;

% 2 x CO2
% Ts = 7.919; %C
% T0 = 3.87; %C
% lambda = 1.62e7; %m^6kg^-1s^-1
% mu = 22e-8; %deg s m^-3

Q = lambda.*(alpha.*(Ts-T0)+beta.*(X(1,:)./100-SS./100))/(1+lambda*alpha*mu);


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


%% Bifurcation analysis- nontrivial eq
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.001,'h_max',0.002,'h_min',0.0001,'ItMX',2000);
bd = coco(prob,'eq1_CO2','ode','isol','ep',... 
    @(x,p)BoxModel_2DH_coco(x,p,'FamousB_2xCO2'),xeq,{'H' 'gamma'},[H gamma],...  
    'H',[-0.6 0.6]);

% bd = coco_bd_read('eq1');

% 2 x CO2
% bd = coco_bd_read('eq1_CO2');

figure(3); hold on;
[SN1,HB1,BP1] = bifdiag1D(bd,1,'H');
set(gca,'FontSize',16)
xlabel('H (Sv)');
ylabel('S_{N} (scaled)');
box on;

%% Continuing periodic orbits
prob = coco_prob();
prob = coco_set(prob,'cont','NAdapt',1,'PtMX',50,'ItMX',[4000 0],'h0',0.01,'h_max',0.5,'h_min',0.01);
bd_po = coco(prob,'po1_CO2','ode','HB','po',... %branch, toolbox family, initial point, point type/branch type
    'eq1_CO2', HB1,... primary branch name, solution label
    {'H'},{[-0.6 0.6]}); %continuation parameter, continuation domain

% bd_po = coco_bd_read('po1');

% 2 x CO2 
% bd_po = coco_bd_read('po1_CO2');

figure(3); hold on;
FP1 = plotPO1D(bd_po,1,'H');

%% Zoom for plot
% prob = coco_prob();
% prob = coco_set(prob,'cont','h0',1e-6,'h_max',1e-3,'h_min',1e-6,'ItMX',2000);
% bd2 = coco(prob,'eq1_zoom','ode','ep','ep',... 
%     'eq1',48,...  
%     'H',[0.21 0.4]);

% bd2 = coco_bd_read('eq1_zoom');

figure(4); hold on;
[SN2,HB2,BP2] = bifdiag1D(bd2,1,'H');

prob = coco_prob();
prob = coco_set(prob,'cont','NAdapt',1,'PtMX',50,'ItMX',[200 0],'h0',1e-8,'h_max',1e-8,'h_min',1e-8);
bd_po2 = coco(prob,'po1_zoom','ode','HB','po',... %branch, toolbox family, initial point, point type/branch type
    'eq1_zoom', HB2,... primary branch name, solution label
    {'H'},{[0.2 0.4]}); %continuation parameter, continuation domain

% bd_po2 = coco_bd_read('po1_zoom');

figure(4); hold on;
plotPO1D(bd_po,1,'H');
box on;
set(gca,'FontSize',16)
xlim([0.2 0.23])

%% 2 Parameter Continuation

%% Continuing Hopf point
% prob = coco_prob();
% prob = coco_set(prob,'cont','h0',0.001,'h_max',0.005,'h_min',0.0001,'ItMX',400);
% bd_hb = coco(prob,'hb1','ode','HB','HB',... 
%    'eq1',HB1,...
%    {'H' 'gamma'},{[-0.6 0.6] [0 2]});

bd_hb = coco_bd_read('hb1');

% 2 x CO2
% bd_hb = coco_bd_read('hb1_CO2');

%% Continuing Fold points
% bd_sn1 = coco(prob,'sn1','ode','SN','SN',... 
%    'eq1',SN1(1),...  
%    {'H' 'gamma'},{[-0.6 0.6] [0 2]});
% 
% bd_sn2 = coco(prob,'sn2','ode','SN','SN',... 
%    'eq1',SN1(2),...  
%    {'H' 'gamma'},{[-0.6 0.6] [0 2]});

bd_sn1 = coco_bd_read('sn1');
bd_sn2 = coco_bd_read('sn2');

% 2 x CO2
% bd_sn1 = coco_bd_read('sn1_CO2');
% bd_sn2 = coco_bd_read('sn2_CO2');

%% 2 Param bifurcation diagram

br = [160 82 45]./255;
dg = [77 149 66]./225;

H_hb = coco_bd_col(bd_hb,'H');
G_hb = coco_bd_col(bd_hb,'gamma');

H_sn1 = coco_bd_col(bd_sn1,'H');
G_sn1 = coco_bd_col(bd_sn1,'gamma');

H_sn2 = coco_bd_col(bd_sn2,'H');
G_sn2 = coco_bd_col(bd_sn2,'gamma');

BT = [0.2268,0.1559];

figure;
plot(H_hb,G_hb,'color',br,'Linewidth',3);
hold on
plot(H_sn1,G_sn1,'color',dg,'Linewidth',3);
plot(H_sn2,G_sn2,'color',dg,'Linewidth',3);
plot(BT(1),BT(2),'k*','Linewidth',1,'Markersize',8);
xlabel('H')
ylabel('\gamma')
axis([0.15 0.25 0 1])
set(gca,'fontsize',16)

%% Zoom for plot
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.001,'h_max',0.0015,'h_min',0.0005,'ItMX',2000);
bd_z = coco(prob,'eq1_zoom','ode','isol','ep',... 
    @BoxModel_2DH_coco,xeq,{'H' 'gamma'},[H gamma],...  
    'H',[0 0.15]);

%% Plotting

H_1 = coco_bd_col(bd,'H');
X_1 = coco_bd_col(bd,'x');

Q_1 = lambda.*(alpha.*(Ts-TN)+beta.*((X_1(1,:)./100+S0)-SS));

figure(5); clf; hold on;
plot(H_1,Q_1,'b','Linewidth',3)
set(gca,'FontSize',16)
title('Bifurcation diagram: H vs Q');
xlabel('H (Sv)');
ylabel('Q (Sv)');

figure;
hold on;
bifdiag1D(bd_z,1,'H');
plotPO1D(bd_po,1,'H');

figure;
subplot(1,2,1); hold on;
bifdiag1D(bd,1,'H');
plotPO1D(bd_po,1,'H');
box on;
set(gca,'FontSize',16)
xlabel('H')
ylabel('S_{N} (scaled)')
subplot(1,2,2); hold on;
bifdiag1D(bd,2,'H');
plotPO1D(bd_po,2,'H');
box on;
set(gca,'FontSize',16)
xlabel('H')
ylabel('S_{T} (scaled)')