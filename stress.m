function [Gstress,Lstress] = stress(theta,Q_bar,strain_g)
c=cosd(theta);
s=sind(theta);
T=[c^2 s^2 2*s*c;
   s^2 c^2 -2*s*c;
   -s*c c*s c^2-s^2 ];
Gstress = Q_bar*strain_g;
Lstress = T*Gstress;
end