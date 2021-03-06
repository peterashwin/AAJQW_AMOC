function [ z ] = BoxModel_2DH_coco( x,p,GCM )
% 3-box function for COCO bifurcation software 
%  GCM - string stating which GCM parameters to use (see below)


load('parameters.txt');

if strcmp(GCM,'FamousA') == 1
    pidx = 1;
elseif strcmp(GCM,'FamousB_1xCO2') == 1
    pidx = 2;
elseif strcmp(GCM,'FamousB_2xCO2') == 1
    pidx = 3;
elseif strcmp(GCM,'HadGEM2_1xCO2') == 1
    pidx = 4; 
elseif strcmp(GCM,'HadGEM2_2xCO2') == 1
    pidx = 5;
elseif strcmp(GCM,'HadGEM2_4xCO2') == 1
    pidx = 6;  
else
    fprintf('\n error : invalid GCM string\n\n')
    return
end

S0 = 0.035; % dimensionless

H = p(1,:);

alpha = 0.12; %kg m^-3 C^-1
beta = 790.0;%kg m^-3 
Y = 100*3.15e7; %sec/year

SN_eq = 0.034912;
ST_eq = 0.035435;
SS_eq = 0.034427;
SIP_eq = 0.034668;
SB_eq = 0.034538;

VN = parameters(1,pidx)*1e16;%m^3
VT = parameters(2,pidx)*1e16; %m^3
VS = parameters(3,pidx)*1e16; %m^3
VIP = parameters(4,pidx)*1e16; %m^3
VB = parameters(5,pidx)*1e16; %m^3

Ts = parameters(14,pidx); %C
T0 = parameters(15,pidx); %C

C = VN*SN_eq+VT*ST_eq+VS*SS_eq+VIP*SIP_eq+VB*SB_eq; %m^3 (total sum of V*S_eq)

lambda = parameters(17,pidx)*1e7; %m^6kg^-1s^-1
%gamma = parameters(22,pidx);
gamma = p(2,:);
mu = parameters(16,pidx)*1e-8; %deg s m^-3

KN = parameters(18,pidx)*1e6; %m^3s^-1
KS = parameters(19,pidx)*1e6; %m^3s^-1

FN = (parameters(10,pidx)+parameters(6,pidx).*H).*1e6;
FT = (parameters(12,pidx)+parameters(7,pidx).*H).*1e6; %m^3s^-1


SN = x(1,:);
ST = x(2,:);
SS = (0.034427-S0).*100;
%SIP = (0.034668-S0).*100;
SB = (0.034538-S0).*100;
SIP = 100.*(C-(VN.*SN+VT.*ST+VS.*SS+VB.*SB)./100-S0.*(VB+VN+VT+VIP+VS))./VIP;

q = lambda.*(alpha.*(Ts-T0)+beta.*(x(1,:)./100-SS./100))/(1+lambda*alpha*mu);
aq = abs(q);


z1p = (Y./VN).*(q.*(ST./100-SN./100)+KN.*(ST./100-SN./100)-FN.*S0);
   
z2p = (Y./VT).*(q.*(gamma.*SS./100+(1-gamma).*SIP./100-ST./100)+...
    KS.*(SS./100-ST./100)+KN.*(SN./100-ST./100)-FT.*S0);
    

z1n = (Y./VN).*(aq.*(SB./100-SN./100)+KN.*(ST./100-SN./100)-FN.*S0);
    
z2n = (Y./VT).*(aq.*(SN./100-ST./100)+KS.*(SS./100-ST./100)...
    +KN.*(SN./100-ST./100)-FT.*S0);


z(1,:) = z1p.*(q>=0)+z1n.*(q<0);
z(2,:) = z2p.*(q>=0)+z2n.*(q<0);

end