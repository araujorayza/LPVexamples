
%% Config iniciais
clear all, clc; close all;

%% Tolerancia de converg�ncia
tol = -1e-9;

%% Main
s = tf('s');
cont = 'PI';

%% Sistema 1

% A = [0 1 ; 0 -1 ];
% B = [0 1]';
% D = 0;
% J = [1 0 ;0 1 ];

% Gn = 1/(s+1);


%% Sistema POLITOPICO
clc
clear all

thetamin=0;
thetamax=1;

A0=[0 1 0;0 0 1;0  0 -1];
A1=[0 0 0;0 0 0;0 -1  0];
A2=[0 0 0;0 0 0;0  0 -1];
A3=[0 0 1;0 0 0;0  0 -1];


Bp = [0 0 1]';
Cp = [1 0]; 
Dp = 0;

Jp = [1 0;1 1; 0 1];

% VARIAVEIS DO ROLMIP
Ap{1} = {1,A0};
Ap{2} = {1,A1};
bounds=[thetamin thetamax;
        thetamin thetamax];
label='A';

A=rolmipvar(Ap,label,bounds)
B=rolmipvar(Bp,'B');
C=rolmipvar(Cp,'C');
D=rolmipvar(Dp,'D');
J=rolmipvar(Jp,'J');
% Sistema 3


[n,m] = size(B);
degP=1;
% Vari�veis LMI
%S = rolmipvar(n, n ,'S','symmetric',{[1 0],[0 1]},bounds);
S=rolmipvar(n,n,'S','sym',[2],[degP]);
Y = rolmipvar(m, n ,'Y','full',{[1 0],[0 1]},bounds);
Ro=rolmipvar('Ro','scalar');
eye=rolmipvar(eye(m),'eye');

T11 = S*A' + A*S - Y'*B' - B*Y ;
T12 = B;
T13 = Y';

T21 = B';
T22 = (-1)*Ro*eye;
T23 = zeros(m)*Ro;

T31 = Y;
T32 = zeros(m)*Ro;
T33 = (-1)*eye;

T = [T11 T12 T13;
     T21 T22 T23;
     T31 T32 T33];
LMIs = [S > 0];
LMIs = [LMIs, T < 0];

%L = [S > 0, T < 0 ];

optimize(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'))

%%
sol = solvesdp(LMIs,Ro)
[p,d] = checkset(LMIs)
clc
if(min(p)>tol)
    polS = double(S);
    polY = double(Y);
    %K = ((inv(Ss)*Ys')'/J')';
    K = ((inv(polS)*polY')'/J')';
     double(Ro)
else
     disp('infactivel')     
end


%%

c = @(x) [real(eig(0.05*eye(n)+A-B*(J*[x(1);x(2)])'))]; % restri��es inequa��es 
ceq = @(x)[];

nonlinfcn = @(x)deal(c(x),ceq(x));
obj = @(x) mse(J*[x(1);x(2)],inv(polS)*polY'); % obj: erro m�dio quadr�tico ser m�nimo
opts = optimoptions(@fmincon,'Algorithm','sqp');
z1 = fmincon(obj,[0;0],[],[],[],[],[],[],nonlinfcn,opts);
[cout,ceqout] = nonlinfcn(z1) % verifica se satisfez as inequa��es 


%%

Kvelho = K; 
K = z1;



%%

Ab = A-B*(J*K)';
Bb = B;
Cb = (J*K)'; 
Db = 0; 


MF_sol = Cb*inv(s*eye(n)-Ab)*Bb+Db 

Cpid = K(2)+K(1)/s;


%% Simulink

%%

MF_real = feedback(Cpid*Gn,1)


impulse(MF_sol)
hold on 
impulse(MF_real)

figure 

step(MF_sol)
hold on 
step(MF_real)

%% 
%clc; close all;

%% Display
double(Ro)

