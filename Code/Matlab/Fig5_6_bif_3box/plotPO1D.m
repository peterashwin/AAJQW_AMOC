function [ FP_lab ] = plotPO1D( bd,var,par )
%Plots 1D bifurcation diagram - COCO
%   bd - branch run ; var - variable index to plot ; par - cont par name

hold on

dg = [77 149 66]./225;

% get min and max of PO
Xmax = coco_bd_col(bd,'MAX(x)');
Xmin = coco_bd_col(bd,'MIN(x)');
p = coco_bd_col(bd,par);


%stability
stab = coco_bd_col(bd,'po.test.USTAB');
changestab = [1];

for i = 2:length(stab)
    if stab(i) ~= stab(i-1)
        changestab = [changestab,i];
    end
end

for j = 1:length(changestab)-1
    ind_change = changestab(j);
    if stab(ind_change) == 0
        plot(p(ind_change:changestab(j+1)),Xmax(var,ind_change:changestab(j+1)),'m',...
            'LineWidth',2);
        plot(p(ind_change:changestab(j+1)),Xmin(var,ind_change:changestab(j+1)),'m',...
            'LineWidth',2)
    else
        plot(p(ind_change:changestab(j+1)),Xmax(var,ind_change:changestab(j+1)),'m--',...
            'LineWidth',2);
        plot(p(ind_change:changestab(j+1)),Xmin(var,ind_change:changestab(j+1)),'m--',...
            'LineWidth',2)
    end
end
       
ind_change = changestab(end);
if stab(ind_change) == 0
    plot(p(ind_change:end),Xmax(var,ind_change:end),'m','LineWidth',2);
    plot(p(ind_change:end),Xmin(var,ind_change:end),'m','LineWidth',2)
else
    plot(p(ind_change:end),Xmax(var,ind_change:end),'m--','LineWidth',2);
    plot(p(ind_change:end),Xmin(var,ind_change:end),'m--','LineWidth',2)
end

% Plot special points using indices
%extract labels
SN = coco_bd_idxs(bd,'SN');

FP_lab = coco_bd_labs(bd,'SN');

plot(p(SN),Xmax(var,SN),'o','color',dg,'LineWidth',2,'MarkerSize',8);
plot(p(SN),Xmin(var,SN),'o','color',dg,'LineWidth',2,'MarkerSize',8);

end

