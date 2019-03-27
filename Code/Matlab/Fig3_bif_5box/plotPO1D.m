function [ FP_lab ] = plotPO1D( bd,var,par,run )
%Plots 1D bifurcation diagram - COCO
%   bd - branch run ; var - variable index to plot ; par - cont par index!!
%   run - run name
hold on

dg = [77 149 66]./225;

LABS = coco_bd_labs(bd);

% get min and max of PO
Xmin = NaN(1,length(LABS));
Xmax = NaN(1,length(LABS));
p = NaN(1,length(LABS));

for i = 1:length(LABS)
    sol = po_read_solution(run,i);
    Xmin(LABS(i)) = min(sol.xbp(:,var),[],1);
    Xmax(LABS(i)) = max(sol.xbp(:,var),[],1);
    p(LABS(i)) = sol.p(par);
end

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
        plot(p(ind_change:changestab(j+1)),Xmax(ind_change:changestab(j+1)),'m',...
            'LineWidth',2);
        plot(p(ind_change:changestab(j+1)),Xmin(ind_change:changestab(j+1)),'m',...
            'LineWidth',2)
    else
        plot(p(ind_change:changestab(j+1)),Xmax(ind_change:changestab(j+1)),'m--',...
            'LineWidth',2);
        plot(p(ind_change:changestab(j+1)),Xmin(ind_change:changestab(j+1)),'m--',...
            'LineWidth',2)
    end
end
       
ind_change = changestab(end);
if stab(ind_change) == 0
    plot(p(ind_change:end),Xmax(ind_change:end),'m','LineWidth',2);
    plot(p(ind_change:end),Xmin(ind_change:end),'m','LineWidth',2)
else
    plot(p(ind_change:end),Xmax(ind_change:end),'m--','LineWidth',2);
    plot(p(ind_change:end),Xmin(ind_change:end),'m--','LineWidth',2)
end

% Plot special points using indices
%extract labels
SN = coco_bd_idxs(bd,'SN');

FP_lab = coco_bd_labs(bd,'SN');

plot(p(SN),Xmax(:,SN),'o','color',dg,'LineWidth',2,'MarkerSize',8);
plot(p(SN),Xmin(:,SN),'o','color',dg,'LineWidth',2,'MarkerSize',8);

end

