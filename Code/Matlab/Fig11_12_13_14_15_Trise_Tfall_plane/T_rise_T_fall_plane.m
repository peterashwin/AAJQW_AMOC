%% (T_rise, T_fall) - plane: the scan

clear

H0    = 0;
Hpert = 0.5; %0.37; 
T0    = 0;
Tpert = 150;

bd = coco_bd_read('eq1_CO2');  %bifurcation data 

H = coco_bd_col(bd,'H');
X = coco_bd_col(bd,'x');
SN = coco_bd_idxs(bd,'SN');
HB = coco_bd_idxs(bd,'HB');

Son1  = interp1(H(1:SN(1)),X(1,1:SN(1)),H0,'spline');
Son2  = interp1(H(1:SN(1)),X(2,1:SN(1)),H0,'spline');
Soff1 = interp1(H(SN(2):end),X(1,SN(2):end),H0,'spline');
Soff2 = interp1(H(SN(2):end),X(2,SN(2):end),H0,'spline');


Son        = [Son1; Son2];
Soff       = [Soff1; Soff2];
Tfall      = 0:10:200;
Trise      = 0:10:200;
AA         = NaN(length(Tfall), length(Trise));

for n = 1 : length(Tfall)
    TFALL = Tfall(n);
    for m = 1: length(Trise)
        [ AA(n,m),~,~] = Rtippingindicator(...
            H0,Hpert,T0,Tpert,Trise(m),TFALL,Son,Soff);
    end
    disp (n)
end
% The boundary curve


for n = 1:length(Trise)
    testfun =@(TT)(Rtippingindicator(H0,Hpert,T0,Tpert,...
        Trise(n),TT,Son,Soff));
    a       = Trise(1);
    b       = Trise(end);
    diff    = 1;
    if testfun(b) == -1
        Tfall1(n) = Trise(end)+10;
    elseif testfun(a) == 1
        Tfall1(n) = Trise(1)-10;
    else
        while(diff >= 1e-1)
            p=(a+b)/2;
            fp = testfun(p);
            fa = testfun(a);
            fb = testfun(b);
            if fa*fp<0
                b = p;
            else
                a = p;
            end
            diff = b-a;
        end
        Tfall1(n)=a;
    end
    disp(n)
end
% Trise = (a+b)/2;
% figure
% plot(Trise,Tfall1)
% axis([0 600 0 600])

%
h = figure;
set(0,'defaulttextInterpreter','latex') %latex axis labels
contourf(Trise,Tfall,AA,'LineStyl','none')
hold on
plot(Trise,Tfall1,'r-','LineWidth',2)
set(gca,'FontSize',18)
xlabel('$T_{rise}$')
ylabel('$T_{fall}$')
title(['$T_{pert} =$' num2str(Tpert)])
colormap([1 1 1; .7 .7 .7])
axis([Trise(1) Trise(end) Tfall(1) Tfall(end)])
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
    [pos(3), pos(4)])
print(h,'TrisVSTfall_Btipping_Hpert05_Tpert150','-dpng','-r0')
print(h,'TrisVSTfall_Btipping_Hpert05_Tpert150','-dpdf','-r0')