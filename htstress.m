function [Rstress,Lhtstress] = htstress(theta,dt,dc,Q_bar,strain_ght,alpha_g,beta_g)
c=cosd(theta);
s=sind(theta);
T=[c^2 s^2 2*s*c;
   s^2 c^2 -2*s*c;
   -s*c c*s c^2-s^2 ];
tstrain = dt*alpha_g;
hstrain = dc*beta_g;
rstrain = strain_ght - tstrain- hstrain;
Rstress= Q_bar*rstrain;
Lhtstress = T*Rstress;
end