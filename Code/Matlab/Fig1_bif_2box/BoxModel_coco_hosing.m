function [ z ] = BoxModel_coco_hosing( x,p )

eta1 = 3;
eta3 = 0.3;

eta2 = p;

S = x(1,:);
T = x(2,:);




z1p = eta2 - S.*( eta3 + abs(T-S) );
z2p = eta1 - T.*(  1 + abs(T-S) );
    

z(1,:) = z1p;
z(2,:) = z2p;

end