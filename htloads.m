function [NC,MC] = htloads(dc,Q_bar,h2,h1,beta_g,NC,MC) 
NC = NC + dc*Q_bar*beta_g*(h2-h1);
MC = MC + (dc/2)*Q_bar*beta_g*(h2^2-h1^2);
end