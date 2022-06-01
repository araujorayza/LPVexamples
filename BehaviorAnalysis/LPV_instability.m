clc;
close all;
clear all;

% The system here is 
% xdot = [ 0                 1;
%         -1 - 0.5*rho(t) -0.2]
syms t rho

% Note that the eigenvalues of the matrix A are negative for fixed values
% of rho(t) between -1 and 1

eigenvalues = eig([           0,   1;
                    -1 - 0.5*rho, -0.2])

eig1 = eigenvalues(1);
eig2 = eigenvalues(2);

hold on
rho_interval = [-1 1];
fplot(real(eig1),imag(eig1),rho_interval)
fplot(real(eig2),imag(eig2),rho_interval)
xlabel('real part of eigenvalues')
ylabel('imaginary part of eigenvalues')
zlabel('rho')
legend('eigenvalue1','eigenvalue2')
grid on
%% to be continued

% But if we choose a time varying rho(t) between -1 and 1, 
% say, rho(t) = cos(2t)
time=0:0.1:7;

dx = @(t,x) [               0,   1;
            -1 - 0.5*cos(2*t), -0.2]*x;
          
rng('shuffle');   %seeds the random number generator based on the current
%     time so that RAND, RANDI, and RANDN produce a different sequence of
%     numbers after each time you call rng    

xin = -1 + 2*rand([2 7]) %setup 7 different initial conditions

Y = cell(length(xin));  %we will save the trajectory of te system for each
                        % initial condition in a cell
    
for i=1:length(xin)
    x0=xin(:,i);    %for each initial condition
    [t,y] = ode45(dx,time,x0);  %solve the ODE for that initial condition
    Y{i}=y; %save the result in the appropriate cell
end 

h=figure; hold on;
for i=1:length(xin) %plot them all in the same graph
    plot(xin(1,i),xin(2,i),'*',Y{i}(:,1),Y{i}(:,2));
    legend('Initial Position','Trajectory')
end 
grid on;
print(h,'LPV_instability','-dpdf')