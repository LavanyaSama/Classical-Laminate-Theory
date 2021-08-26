function [A,B,D] = abd(Q_bar,h2,h1,A,B,D) 
A = A + Q_bar*(h2-h1);
B = B + 0.5*(Q_bar*(h2^2-h1^2));
D = D +(1/3)*(Q_bar*(h2^3-h1^3));
end