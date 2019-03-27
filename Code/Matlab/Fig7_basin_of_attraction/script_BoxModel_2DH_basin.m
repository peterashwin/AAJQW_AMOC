%% Richard Wood box model - Basins of attraction 

clear;
close all;
%%  interpolating the branches

bd = coco_bd_read('eq1_CO2'); % bifurcation data eq1_CO2 for 2*CO2 and eq1 for 1*CO2


H = coco_bd_col(bd,'H');
X = coco_bd_col(bd,'x');

SN = coco_bd_idxs(bd,'SN');
HB = coco_bd_idxs(bd,'HB');
H_sn1 = H(SN(1));
H_sn2 = H(SN(2));
H_hb  = H(HB);

up_bran1  = @(input)interp1(H(1:SN(1)),X(1,1:SN(1)),...
    input,'spline');
lo_bran1 = @(input)interp1(H(SN(2):end),X(1,SN(2):end),...
    input,'spline');
un_bran1 = @(input)interp1(H(SN(1):SN(2)),X(1,SN(1):SN(2)),...
    input,'spline');

up_bran2  = @(input)interp1(H(1:SN(1)),X(2,1:SN(1)),...
    input,'spline');
lo_bran2 = @(input)interp1(H(SN(2):end),X(2,SN(2):end),...
    input,'spline');
un_bran2 = @(input)interp1(H(SN(1):SN(2)),X(2,SN(1):SN(2)),...
    input,'spline');
%% Calculating the basins

stepsize=10;

Sn_min = -0.15;
Sn_max =  0.05;

St_min = -0.05;
St_max =  0.2;

N = 25;
M = 25; % N*M = number of boxes

Sn_basin = linspace(Sn_min,Sn_max,N);
St_basin = linspace(St_min,St_max,M);
tspan    = [0,4000];

% Xends = NaN(N,M,2);
K     = 200;
B     = zeros(K,N,M);
Hbas  = linspace(H(SN(1)),H(SN(2)),K);
Vup   = zeros(1,K);
Vlo   = zeros(1,K);


for k = 1:K
    
    Xu =[up_bran1(Hbas(k));up_bran2(Hbas(k))];
    Xl =[lo_bran1(Hbas(k));lo_bran2(Hbas(k))];
    Xs =[un_bran1(Hbas(k));un_bran2(Hbas(k))];
    
    for m=1:N
        for n=1:M
            x0          = [Sn_basin(n);St_basin(m)];
            [tend,xend] = MyIVP_basin(@(t,x)BoxModel_2DH_IVP(t,x,Hbas(k)),...
                x0,tspan,stepsize,Xu,Xl);
            B(k,m,n)      = tend;
            
            dist1 = norm(xend - Xu);
            dist2 = norm(xend - Xl);
            if dist1 <= 1e-2
                Vup(k) = Vup(k)+1;
                     B(k,m,n)      = tend;
            elseif dist2<=1e-2
                Vlo(k) = Vlo(k)+1;
                     B(k,m,n)      = -tend;
            end
            disp([k,n,m,min(dist1,dist2)])
        end
        
    end
 end
figure;
plot(Hbas,Vlo/(N*M))
hold on;
plot(Hbas,Vup/(N*M))
%% Movie for the basins vs H
V = VideoWriter('MOC_basin_2co2','MPEG-4');
V.FrameRate = 8;
V.Quality = 100;
open(V);
warning off
for k = 80:-1:1
    Xu =[up_bran1(Hbas(k));up_bran2(Hbas(k))];
    Xl =[lo_bran1(Hbas(k));lo_bran2(Hbas(k))];
    Xs =[un_bran1(Hbas(k));un_bran2(Hbas(k))];
    A=B(k,:,:);
    A=reshape(A,[],size(A,2),1);
    A=abs(A);
%     figure;
    subplot(8,1,1:5)
    set(0,'defaulttextInterpreter','latex') %latex axis labels
    hold on
    BASIN = pcolor(Sn_basin,St_basin,A);
    BASIN.LineStyle = 'none';
    colormap hot
    caxis([0 2000]);
    colorbar
    plot(Xu(1),Xu(2),'bo','LineWidth',3)
    plot(Xl(1),Xl(2),'bo','LineWidth',3)
    
    xlim([ Sn_min Sn_max])
    ylim([ St_min St_max])
    
    [t,xt] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,Hbas(k)),...
        [10000 0],Xs-[0.01; 0]);
    plot(xt(13:end,1),xt(13:end,2),'-','LineWidth',3,'Color',[0.5 0.5 0.5])
    [t,yt] = ode45 (@(t,x)BoxModel_2DH_IVP(t,x,Hbas(k)),...
        [10000 0],Xs-[0; 0.0001]);
    plot(yt(:,1),yt(:,2),'-k','LineWidth',3,'Color',[0.5 0.5 0.5])
    plot(Xs(1),Xs(2),'g+','LineWidth',3,'MarkerSize',10)
    xlabel('$S_N$');
    ylabel('$S_T$','Rotation',0);
    set(gca,'FontSize',16)  
    title(['H = ' num2str(Hbas(k))])
    hold off
%     
    subplot(8,1,7:8)
    plot(Hbas,Vup./(N*M),'-b',Hbas,Vlo./(N*M),'-r','LineWidth',3);
    hold on
    plot(Hbas(k)*ones(1,10),linspace(-0.1,1.1,10),'-.k',...
        'LineWidth',2)
    xlim([Hbas(80),Hbas(1)])
    ylim([-0.1 1.1])
    set(gca,'FontSize',16)
    set(gcf,'Position',[0.1e3 0.07e3 .8e3 0.7e3])
    xlabel('H')
    ylabel('Basins Volume')
    hold off

    F(k) = getframe(gcf);
    writeVideo(V,F(k));

 end
close(V);
winopen('MOC_basin_2co2.mp4')


