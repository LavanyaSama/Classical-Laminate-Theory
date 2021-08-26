clear all
clc
%% properties of laminate
E1=181e9;
E2=10.3e9;
v1=0.28;
G12=7.17e9;
v2=(E2/E1)*v1;

%% Elements of Reduced stiffness matrix
Q11=E1/(1-(v1*v2));
Q22=E2/(1-(v1*v2));
Q12=(v1*E2)/(1-(v1*v2));
Q66=G12;
Q16=0;
Q26=0;
%% Load per unit length
N=[1 0 0];% in the order Nx, Ny,Nxy(N/m)
%% Moment per unit length
M=[0 0 0];% in the order Mx, My,Mxy
%% distance from mid axis
h=[-0.0075 -0.0025  0.0025 0.0075];
%% coefficient of thermal expansion
alpha=[2e-8 22.5e-6 0];
%% coefficient of moisture expansion
beta=[0 0.6 0];
Q=[Q11 Q12 0;
   Q12 Q22 0;
   0    0  Q66];
%% Laminate Staking sequence
theta=[0 90 0];% Enter a symmetric matrix
dt=-75;
dc=0;
fid = fopen('Report.txt','w');
fprintf(fid,'Input:\n\n'); 
fprintf(fid,'Laminate Properties-Units(Pa):'); 
fprintf(fid,'\n\n');
fprintf(fid,'E1 = \t');
fprintf(fid,'%d\n',E1); 
fprintf(fid,'E2 = \t');
fprintf(fid,'%d\n',E2); 
fprintf(fid,'v1 = \t');
fprintf(fid,'%.3f\n',v1); 
fprintf(fid,'G12 = \t');
fprintf(fid,'%d\n',G12); 
fprintf(fid,'Staking Sequence(degrees) = \t'); 
fprintf(fid,'%i\t',theta);
fprintf(fid,'\n\n');
fprintf(fid,'Load per unit length(Nx,Ny,Nxy) = \t'); 
fprintf(fid,'%i\t',N);
fprintf(fid,'\n\n');
fprintf(fid,'Coefficients of thermal expansion(alpha x,alpha y,alpha xy) = \n'); 
fprintf(fid,'%f\t',alpha);
fprintf(fid,'\n\n');
fprintf(fid,'Coefficients of moisture expansion(beta x,beta y,beta xy) = \n'); 
fprintf(fid,'%f\t',beta);
fprintf(fid,'Input:\n\n'); 
fprintf(fid,'Reduced Stiffness matrix:\n');
for i=1:3
for j=1:3
fprintf(fid,'%e\t',Q(i,j));
end
fprintf(fid,'\n');
end
fprintf(fid,'\n\n');
%% Loop for calculating  Reduced transformed stiffness  matrix, coefficients of thermal and moisture expansion in global axes
for i=1:length(theta)  
c=cosd(theta(i));
s=sind(theta(i));
fprintf(fid,'Angle = \t');  
fprintf(fid,'%d\n\n',theta(i));
T=[c^2 s^2 2*s*c;
    s^2 c^2 -2*s*c;
    -s*c c*s c^2-s^2 ];
R=[1 0 0;
    0 1 0;
    0 0 2];
% Reduced transformed stiffness  matrix
fprintf(fid,'Transformed Reduced Stiffness matrix:\n\n');
Q_bar(:,:,i)=inv(T)*Q*R*T*inv(R);
for m=1:3
for n=1:3
fprintf(fid,'%e\t',Q_bar(m,n,i));
end
fprintf(fid,'\n');
end
% alpha in global axis
alpha_g(:,:,i) = R*inv(T)*alpha';
% beta in global axis
beta_g(:,:,i) = R*inv(T)*beta';
end
for j=1:length(theta)
 %% Initialisation of A,B and D matrix
A=zeros(1);
B=zeros(1);
D=zeros(1);
NT=zeros(3,1);
NC=NT;
MT=zeros(3,1);
MC=MT;
SR=zeros(1);
SR1=zeros(1);
%% Calculating A,B and D matrices and hygrothermal loads and matrices
for i=1:length(theta)
[A,B,D] = abd(Q_bar(:,:,i),h(i+1),h(i),A,B,D);  
[NT,MT] = tloads(dt,Q_bar(:,:,i),h(i+1),h(i),alpha_g(:,:,i),NT,MT);
[NC,MC] = htloads(dc,Q_bar(:,:,i),h(i+1),h(i),beta_g(:,:,i),NC,MC);
end
%% Printing A,B and D matrices
fprintf(fid,'\n');
fprintf(fid,'A matrix:\n');
for m=1:3
for n=1:3
fprintf(fid,'%e\t',A(m,n));
end
fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'B matrix:\n');
for m=1:3
for n=1:3
fprintf(fid,'%e\t',A(m,n));
end
fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'D matrix:\n');
for m=1:3
for n=1:3
fprintf(fid,'%e\t',A(m,n));
end
fprintf(fid,'\n');
end
fprintf(fid,'\n');
%% Printing Thermal loads
fprintf(fid,'Thermal load per unit length:\n');
for m=1:3
fprintf(fid,'%e\n',NT(m));
end
fprintf(fid,'Thermal Moment per unit length:\n');
for m=1:3
fprintf(fid,'%e\n',MT(m));
end
%% Terminating loop if all plies damaged
if A==0
    break;
end
%% Combined Stiffness matrix
Matrix = [A,B;
          B,D];
TLoad = [NT;MT];% Thermal
MLoad = [NC;MC];% Hygro
ELoad = [N M];% External load
Load = TLoad + MLoad;
%% Calculating mid plane strains and Plate curvatures
midstrain = inv(Matrix)*ELoad'; 
midstrainht = inv(Matrix)*Load;
% Midplane strains and curvatures
[m1,cu]= straincurv(midstrain);
[mht,cuht]=straincurv(midstrainht);
% distance of each lamina from mid axis
z=[-0.0075 -0.0025  0.0025];
%% Printing miplane strains and plate curvatures
fprintf(fid,'midplane strains and curvatures:\n');
for m=1:6
fprintf(fid,'%e\n',midstrain(m));
end
fprintf(fid,'midplane strains and curvatures for hygrothermal loads:\n');
for m=1:6
fprintf(fid,'%e\n',midstrainht(m));
end
fprintf(fid,'Stresses:\n\n');
%% Calculating Local and Global stresses 
for i=1:length(theta)
    fprintf(fid,'For angle %i :\n',theta(i));
    strain_g(:,:,i) = m1 + z(i)*cu;
    [Gstress,Lstress] = stress(theta(i),Q_bar(:,:,i),strain_g(:,:,i));
    fprintf(fid,'Local Stresses\n\n');
    for m=1:3
    fprintf(fid,'%e\n',Lstress(m));
    end
    fprintf(fid,'Global Stresses\n\n');
    for m=1:3
    fprintf(fid,'%e\n',Gstress(m));
    end
    LS(:,:,i) = Lstress;
    strain_ght(:,:,i) = mht + z(i)*cuht;
    [Rstress,Lhtstress] = htstress(theta(i),dt,dc,Q_bar(:,:,i),strain_ght(:,:,i),alpha_g(:,:,i),beta_g(:,:,i)); 
    fprintf(fid,'Residual local Stresses\n\n');
    for m=1:3
    fprintf(fid,'%e\n',Lhtstress(m));
    end
    fprintf(fid,'Resudual global Stresses\n\n');
    for m=1:3
    fprintf(fid,'%e\n',Rstress(m));
    end
    LS(:,:,i) = Lstress;
    RS(:,:,i) = Lhtstress;
end
%% Failure Analysis
fprintf(fid,'Failure Analysis\n');
% Lamina Srength Properties
ULT = 1500e6;% Ultimate Longitudinal Tensile strength
ULC = 1500e6;% Ultimate Longitudinal Compressive strength
UTT = 40e6;% Ultimate Transverse Tensile strength
UTC = 246e6;% Ultimate Transverse Tensile strength
US=68e6;% Ultimate Shear strength
%% Tsai-Wu Coefficients and strength ratios. Solving quadratic equation obtained from Tsai-wu failure equation "f1.S1+f2.S2+f6.f6+f11.S1^2+f22.S2^2+f66.S6^2+2f12.S1.S2<1" to get strength ratios.
[F] = tsaiwucoef(ULT,ULC,UTT,UTC,US);
    fprintf(fid,'Tsai-wu coefficients\n\n');
    for m=1:6
    fprintf(fid,'%e\t',F(m));
    end
    fprintf(fid,'\n');
for i = 1:length(theta)
    LSS = LS(:,:,i);
    RSS = RS(:,:,i);
    a = F(3)*LSS(1,1)^2 + F(4)*LSS(2,1)^2 + F(5)*LSS(3,1)^2 + 2*F(6)*LSS(1,1)*LSS(2,1);
    b = F(1)*LSS(1,1) + F(2)* LSS(2,1);
    %a1 = F(3)*LSS(1,1)^2 + F(4)*LSS(2,1)^2 + F(5)*LSS(3,1)^2 + 2*F(6)*LSS(1,1)*LSS(2,1);
    b1 = F(1)*LSS(1,1) + F(2)* LSS(2,1) + 2*F(3)*RSS(1,1)*LSS(1,1) + 2*F(4)*RSS(2,1)*LSS(2,1) + 2*F(5)*RSS(3,1)*LSS(3,1) + 2*F(6)*( RSS(1,1)*LSS(2,1) + RSS(2,1)*LSS(1,1) );
    c1 = F(1)*RSS(1,1) + F(2)* RSS(2,1) + F(3)*RSS(1,1)^2 + F(4)*RSS(2,1)^2 + F(5)*RSS(3,1)^2 + 2*F(6)*RSS(1,1)*RSS(2,1);
    if (a~=0 || b~=0)
    syms x
    eq=solve(a*x^2 + b*x - 1); % Solving quadratic equation for obtaining Strength ratios of normal loading
    end 
    if (a~=0 || b1~=0)
    syms x
    eq1=solve(a*x^2 + b1*x + c1 - 1); % Solving quadratic equation for obtaining Strength ratios of Hygrothermal loading
    end
    eq = double(eq);
    eq1 = double (eq1); 
    SR(i) = max(eq);
    SR1(i) = max(eq1);
end
fprintf(fid,'Strength ratios for Normal load\n\n');
for m=1:3
fprintf(fid,'%e\t',SR(m));
end
fprintf(fid,'\n');
fprintf(fid,'Strength ratios for Thermal load\n\n');
for m=1:3
fprintf(fid,'%e\t',SR1(m));
end
fprintf(fid,'\n\n');
%% Ply failure loads
[FPF,index] = min(SR);
[TFPF,index1] = min(SR1);
PF(j)=FPF;
TPF(j)=TFPF;
Angle = theta(index);
for k = 1:length(theta)
if theta(k)==Angle
Q_bar(:,:,k)=zeros(3);% loop to make Q_bar zero to degrade respective damaged ply
end
end
Angle1 = theta(index1);
for l = 1:length(theta)
if theta(l)==Angle1
Q_bar(:,:,l)=zeros(3);% loop to make Q_bar zero to degrade respective damaged ply
end
end
fprintf(fid,'Next step in Complete degradation\n\n\n');
end
fprintf(fid,'First ply failure(N/m) = %f\n',PF(1));
fprintf(fid,'Last ply failure(N/m) = %f\n',PF(end));
fprintf(fid,'First ply failure(N/m) = %f\n',TPF(1));
fprintf(fid,'First ply failure(N/m) = %f\n',TPF(end));
fclose(fid);