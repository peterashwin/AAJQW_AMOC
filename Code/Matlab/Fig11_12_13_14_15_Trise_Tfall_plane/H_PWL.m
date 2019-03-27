function output = H_PWL(t,H0,Hpert,T0,Tpert,Tris,Tfall)

% H0    = 0;
% Hpert = 0.151;
% T0    = 300;
% Tris  = 120;
% Tpert = 100;
% Tfall = 2674;

rris    = (Hpert - H0)/Tris;
rfall   = (H0 - Hpert)/Tfall;

if t <= T0
    output = H0;
elseif and(t>T0,t<=(Tris+T0))
    output = (rris * t + (H0 - (rris * T0))); % linear 
elseif and(t>(Tris+T0), t<=(T0+Tris+Tpert))
    output = Hpert;
elseif and(t>(T0+Tris+Tpert),t<=(T0+Tris+Tpert+Tfall))
    output = rfall * t + (Hpert - (rfall *(T0+Tris+Tpert))); %linear 
elseif t>(T0+Tris+Tpert+Tfall)
    output = H0;
end

end

