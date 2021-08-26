function [NT,MT] = tloads(dt,Q_bar,h2,h1,alpha_g,NT,MT) 
NT = NT + dt*Q_bar*alpha_g*(h2-h1);
MT = MT + (dt/2)*Q_bar*alpha_g*(h2^2-h1^2);
end