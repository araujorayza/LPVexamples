
% Config iniciais
clear all, clc; close all;

% Tolerancia de converg�ncia
tol = -1e-5;

% Main
s = tf('s');

% Sistema

alpha = sdpvar(1);
A = [0 1  0;0 0 1; 0 -1 -2*alpha];
B = [0 0 1]';
C = [1 0]; 
D = 0;
J = [1 0 0;1 1 -1;0 1 1];

%Gn = (s+1)/(s^2+s+1);


% Sistema 3


[n,m] = size(B);

% Vari�veis LMI
S = sdpvar(n, n ,'sym');
Y = sdpvar(m, n ,'full');
sdpvar Ro

T11 = S*A' + A*S - Y'*B' - B*Y ;
T12 = B;
T13 = Y';

T21 = B';
T22 = -Ro*eye(m);
T23 = zeros(m);

T31 = Y;
T32 = zeros(m);
T33 = -eye(m);

T = [T11 T12 T13; T21 T22 T23 ; T31 T32 T33];

LMIs = [S > 0, T < 0, -1.5 <= alpha <= 1.5, uncertain(alpha)];

sol = optimize(LMIs,Ro)
yalmiperror(sol.problem)

[p,d] = checkset(LMIs)
if(min(p)>tol)
    polS = value(S);
    polY = value(Y);
    %K = ((inv(Ss)*Ys')'/J')';
    K = ((inv(polS)*polY')'/J')'
     double(Ro)
else
     disp('infactivel')     
end
%agora como recuperar K?