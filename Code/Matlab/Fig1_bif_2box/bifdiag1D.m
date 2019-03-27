function [ SN_lab,HB_lab,BP_lab ] = bifdiag1D( bd,var,par )
%Plots 1D bifurcation diagram - COCO
%   bd - branch run ; var - variable index to plot ; par - continuation parameter
hold on

br = [160 82 45]./255;
dg = [77 149 66]./225;

%extract columns
p = coco_bd_col(bd,par);
x = coco_bd_col(bd,'x');
xvar = x(var,:);

%stability
stab = coco_bd_col(bd,'ep.test.USTAB');
changestab = [1];

for i = 2:length(stab)
    if stab(i) ~= stab(i-1)
        changestab = [changestab,i];
    end
end

for j = 1:length(changestab)-1
    ind_change = changestab(j);
    if stab(ind_change) == 0
        plot(p(ind_change:changestab(j+1)),xvar(ind_change:changestab(j+1)),'b',...
            'LineWidth',2);
    elseif stab(ind_change) == 1
        plot(p(ind_change:changestab(j+1)),xvar(ind_change:changestab(j+1)),'g',...
            'LineWidth',2);
    else
        plot(p(ind_change:changestab(j+1)),xvar(ind_change:changestab(j+1)),'r',...
            'LineWidth',2);
    end
end
       
ind_change = changestab(end);
if stab(ind_change) == 0
    plot(p(ind_change:end),xvar(ind_change:end),'b','LineWidth',2);
elseif stab(ind_change) == 1
    plot(p(ind_change:end),xvar(ind_change:end),'g','LineWidth',2);
else
    plot(p(ind_change:end),xvar(ind_change:end),'r','LineWidth',2);
end

% Plot special points using indices
%extract labels
SN = coco_bd_idxs(bd,'SN');
HB = coco_bd_idxs(bd,'HB');
BP = coco_bd_idxs(bd,'BP');

SN_lab = coco_bd_labs(bd,'SN');
HB_lab = coco_bd_labs(bd,'HB');
BP_lab = coco_bd_labs(bd,'BP');

plot(p(SN),xvar(:,SN),'o','color',dg,'LineWidth',2,'MarkerSize',8);
plot(p(HB),xvar(:,HB),'d','color',br,'LineWidth',2,'MarkerSize',8);
plot(p(BP),xvar(:,BP),'ko','LineWidth',2,'MarkerSize',8);

end

