function [ indic,rris,rfall ] = Rtippingindicator(H0,Hpert,T0,Tpert,...
    Tris,Tfall,Son,Soff)
% Rtippingindicator: identifies tipping reagions in (Trise Tfall)-plane
%   inputs: the parameter H0,Hpert,T0,Tpert,Tris,Tfall
%   output: -1 if there is no R tipping, if there is R-tipping



rris       = abs(Hpert - H0)/Tris;
rfall      = abs(Hpert - H0)/Tfall;


T          = 1500;%(T0 + Tpert + Tris + Tfall + 1000);
stepsize   = 2;
tspan      = [0 T];

% 
H          = @(tt)H_PWL(tt,H0,Hpert,T0,Tpert,Tris,Tfall);
odefun     = @(tt,x)BoxModel_2DH_IVP(tt,x,H(tt));
[~,~,xend] = MyIVP(odefun,Son,tspan,stepsize);



if norm(xend-Son)<1e-3
    indic = -1;
elseif norm(xend-Soff)<1e-3
    indic = 1;
else
    if norm(xend-Son)<norm(xend-Soff)
        indic = -1;
    else
        indic = 1;
    end
end


end

