% Wrote this script to solve a LPV ODE

clc
clear all

syms t x1(t) x2(t)
assume(t,'real');
a=1.5;

x=[x1;
   x2];

A=[-a*sin(2*t) -a*cos(2*t);
   -a*cos(2*t) a*sin(2*t)];

%can't find an explicit solution
S = dsolve(diff(x) == A*x)

%% Manually finding the state transfer matrix
syms t0 real

X(t)=[1 0;
    0.5*t^2 1];

P(t)=[1 1;
    0.5*t^2 0.5*t^2+2];


simplify(X(t)*inv(X(t0)))
simplify(P(t)*inv(P(t0)))
