% 2 box model 
clear;
close all;

% Initial conditions


S = 1;
T = 1;


eta2 = 0;


% Set up for initial value problem solver for the upper branch

x0 = [S;T];
tspan = [0,1000];
h = 0.1;

[X,t,xeq] = MyIVP(@(t,x)BoxModel_IVP_hosing(t,x,eta2),x0,tspan,h);



% Bifurcation analysis- nontrivial eq 1
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.01,'h_max',0.1,'h_min',0.01,'ItMX',1000);
prob = coco_set(prob,'corr','MX',10);
bd = coco(prob,'eq_2box','ode','isol','ep',... 
    @BoxModel_coco_hosing,xeq,{'eta2'},eta2,...  
    'eta2',[0 2]);
bd   = coco_bd_read('eq_2box');
x    = coco_bd_col(bd,'x');
eta2 = coco_bd_col(bd,'eta2');
SN = coco_bd_idxs(bd,'SN');

psi1   = (  x(2,1:SN)   -  x(1,1:SN)   );
psi2   = (  x(2,SN:end) -  x(1,SN:end) );
psiSN1 = (  x(2,SN)     -  x(1,SN)); 



%
dg = [77 149 66]./225;
figure(3); hold on
plot(eta2(1:SN),psi1,'b',eta2(SN:end),psi2,'g','LineWidth',2)
plot(eta2(SN),psiSN1,'o','color',dg,'LineWidth',2,'MarkerSize',8);

%% ------------------------------------------------------------------------

% the lower brnch 


% Initial conditions


S = 1;
T = 1;


eta2 = 2;


% Set up for initial value problem solver for the upper branch

x0 = [S;T];
tspan = [0,1000];
h = 0.1;

[X,t,xeq] = MyIVP(@(t,x)BoxModel_IVP_hosing(t,x,eta2),x0,tspan,h);



% Bifurcation analysis- nontrivial eq 1
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.01,'h_max',0.1,'h_min',0.01,'ItMX',1000);
prob = coco_set(prob,'corr','MX',10);
bd = coco(prob,'eq_2box','ode','isol','ep',... 
    @BoxModel_coco_hosing,xeq,{'eta2'},eta2,...  
    'eta2',[0 2]);
bd   = coco_bd_read('eq_2box');
y    = coco_bd_col(bd,'x');
eta2 = coco_bd_col(bd,'eta2');


psi3   = (  y(2,:)   -  y(1,:)   );

% psiSN1 = (  x(2,SN)     -  x(1,SN)); 



%
dg = [77 149 66]./225;

plot(eta2,psi3,'b','LineWidth',2)
plot(eta2(1),psi3(1),'o','color',dg,'LineWidth',2,'MarkerSize',8);
