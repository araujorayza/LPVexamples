%Script based on https://yalmip.github.io/example/lpvstatefeedback/
clc
clear all

alpha = sdpvar(1);
Y = sdpvar(3,3);
L = sdpvar(1,3,'full');

A = [0 1 alpha;0 0 1;0 0 0];
B = [0;0;1];
Q = eye(3);
R = 1;


F = [Y >=0];
F = [F, [-A*Y-B*L + (-A*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] > 0]

F = [F, -0.1 <= alpha <= 0.1, uncertain(alpha)];

diag=optimize(F,-trace(Y))

Kpoly = value(L)*inv(value(Y))

cost_poly=trace(value(Y))

%%
clear L F diag

L0 = sdpvar(1,3);
L1 = sdpvar(1,3);
L = L0 + alpha*L1;

F = [Y >=0];
F = [F, [-A*Y-B*L + (-A*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] > 0]

F = [F, -0.1 <= alpha <= 0.1, uncertain(alpha)];

diag=optimize(F,-trace(Y))

%K_GS = (value(L0)+alpha*value(L1))/(value(Y)) %scheduling rule

cost_GS=trace(value(Y))

%%
clear Y F diag

Y0 = sdpvar(3,3);
Y1 = sdpvar(3,3);

% Y1(:,3) = 0;
% Y1(3,:) = 0;
Y = Y0 + alpha*Y1;

F = [Y >=0];
F = [F, [-A*Y-B*L + (-A*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] >= 0]
F = [F, -0.1 <= alpha <= 0.1, uncertain(alpha)];

opt=sdpsettings;
%opt.solver='lmilab';
opt.solver='mosek';

optimize(F,-trace(Y),opt)

cost_Plyap=trace(value(Y)) %this results NaN and it shouldn't
