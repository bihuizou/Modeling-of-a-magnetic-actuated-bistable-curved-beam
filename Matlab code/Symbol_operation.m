clc
clear all
close all

% --- This file is the symbol operation during the modeling process, from
%     which we can obtain the expressions of magnetic torques, energy, and
%     stationary conditions.
%    
% --- Eq1 and Eq2 are the two stationary conditions that need to be solved.
% --- T1expr and T2expr are the expressions of megnatic torques.
%
% --- To save calculation time, we directly substitude the expressions
%     of T1, T2, Eq1 and Eq2 into the multimodal deformation analysis. 

syms a1 a2 x L h pp EE b t II L1 L2 T1 T2 B thetaB m1 m2

y1 = a1*( 1-cos(2*pp*x/L)  );
y2 = a2*( 1-2*x/L - cos( 2.86*pp*x/L ) + 2/2.86/pp*sin(2.86*pp*x/L )  );
y = y1+y2;

dy = diff(y,x);
ddy = diff(y,x,2);

y0 = h/2*( 1-cos(2*pp*x/L)  );
dy0 = diff(y0,x);
ddy0 = diff(y0,x,2); 


ds = int( ( 1/2*(  (dy)^2-(dy0)^2 ) ),x,0,L   );
s0 = int( (1+1/2*(dy0)^2  ),x,0,L   );

p = -EE*b*t*ds/s0;
Uc = -p*ds;
Ub = EE*II/2*int( (  ddy-ddy0  )^2, x,0,L  );

theta1bar = subs(dy0,x,L1  );
theta2bar = subs(dy0,x,L2  );
theta1 = subs(dy,x,L1  );
theta2 = subs(dy,x,L2  );

Um = -T1*(theta1bar-theta1) - T2*(theta2bar-theta2);
U = Um+Uc+Ub;

dUda1 = diff(U,a1);
dUda2 = diff(U,a2);

T1expr = m1*B*sin(theta1-thetaB)
T2expr = m2*B*sin(pp+theta2-thetaB)

Eq1 = subs(subs(dUda1,T1,T1expr),T2,T2expr)
Eq2 = subs(subs(dUda2,T1,T1expr),T2,T2expr)


