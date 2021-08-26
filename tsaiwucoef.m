function [F]=tsaiwucoef(ULT,ULC,UTT,UTC,US)
%Tsai Wu failure criterion:
%h1.S1+h2.S2+h6.S6+h11.S1^2+h22.S2^2+f66.S6^2+2f12.S1.S2

f1=1/(ULT)-1/(ULC);
f2=1/(UTT)-1/(UTC);
f11=1/(ULT*ULC);
f22=1/(UTT*UTC);
f66=1/(US)^2;
f12=-0.5*sqrt(f11*f22);
F=[f1;f2;f11;f22;f66;f12];

end