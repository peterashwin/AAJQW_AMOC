%% Richard Wood box model - variable hosing
clear;
close all;

%% Initial conditions

S0 = 0.035;

% SN = -0.0083;
% ST = 0.0570;

% 2 x CO2
SN = 0.032374702275250;
ST = 0.143450011129485;

% H = 0;


%% Set up for initial value problem solver

x0 = [SN;ST];
tspan = [0,2000];
h = 1;

%% Phase portraits

tstart = 100;
Hshift = 0.5;

Tend = [320:5:345,400,425,450,1400];
%185:1:194;

Colors = [1 0 0;
    1 128/255 0;
    1 230/255 0;
    102/255 204/255 0;
    0 204/255 204/255;
    0 102/255 204/255;
    102/255 0 204/255;
    204/255 0 204/255;
    102/255 51/255 0;
    0 0 0];
set(0,'defaulttextInterpreter','latex') %latex axis labels
figure(1); hold on;

for l= 1:length(Tend)
    tend = Tend(l);
    H=@(t) (t < tstart)*0 + ((tstart<=t)&&(t<=tend))*Hshift + (t>tend)*0;
    [X,t,xeq] = MyIVP(@(t,x)BoxModel_2DH_IVP(t,x,H(t),'FamousB_2xCO2'),x0,tspan,h);
    color = Colors(l,:);
    p(l) = plot(X(1,:),X(2,:),'Linewidth',3,'color',color);
end

leg = repmat({''},length(Tend),1);

for l = 1:length(Tend)
    leg(l) = {num2str(Tend(l)-tstart)};
end

figure(1);
plot(SN,ST,'bo','Linewidth',3,'Markersize',6)
plot(-0.198,0.1546,'bo','Linewidth',3,'Markersize',6)
axis([-0.3 0.1 -0.05 0.4])
legend(p,leg)
set(gca,'FontSize',20)
box on;
xlabel('$S_N$');
ylabel('$S_T$');