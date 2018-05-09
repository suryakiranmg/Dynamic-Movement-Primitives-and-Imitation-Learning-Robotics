clear all; clc; close all;
x(1)=1; alpha_x = 8; beta_z = 6; alpha_z = 25;
c = [ 1 0.6294 0.3962 0.2494 0.1569 0.0988 0.0622 0.0391 0.0246 0.0155];
var=[41.6667 16.5359 6.5359 2.5840 1.0235 0.4054 0.1606 0.0636 0.0252 0.0252]./1000;
w(1)=0;	
load imitation.data
yn = imitation(:,1); y_dn = imitation(:,2); y_ddn = imitation(:,3);

g = ones(1001); g = g(:,1);
f_d = y_ddn - alpha_z*(beta_z*(g - yn) - y_dn);

dt=0.001; t=0:dt:1;
for i = 1:length(t)-1
    x(i+1) = x(i) - alpha_x*x(i)*dt;
end

f_matrix=[];
for i = 1:length(t)
    xdot(i) = - alpha_x * x(i);
    x(i+1) = x(i) + xdot(i) * dt;
    phi=[]; psi=[];
    for j = 1:10
        psi(j)= exp((-1/(2 * var(j))) * (x(i) - c(j))^2 ); end
    for j = 1:10
        phi(j) = (psi(j) * x(i))/sum(psi); end   
    f_matrix = [f_matrix; phi];
end
Computed_weight = ((f_matrix'*f_matrix)^-1)*f_matrix'* f_d

[y,ydot,zdot]=Movt_Planning_Fn(Computed_weight);
figure
subplot(3,1,1); plot(t, yn,t,y); title('Position plot')
subplot(3,1,2); plot(t, y_dn,t,ydot); title('Velocity plot')
subplot(3,1,3);plot(t, y_ddn,t,zdot); title('Acceleration plot')


function [y,ydot,zdot]= Movt_Planning_Fn(Computed_weight)
%Time Constants & Initial Conditions
alpha_z = 25; beta_z = 6; alpha_x = 8;
y0 = 0; x0 = 1; z0 = 0;
N = 10; g=1;
c = [1 0.6294 0.3962 0.2494 0.1569 0.0988 0.0622 0.0391 0.0246 0.0155];
var=[41.6667 16.5359 6.5359 2.5840 1.0235 0.4054 0.1606 0.0636 0.0252 0.0252]./1000;
w = Computed_weight;
dt = 0.001; t = 0:dt:1;
x(1) = x0; xdot(1) = 0; y(1) = y0;
ydot(1) = 0; 
plotvar_f = [];
for i=1:length(t)-1
    for j = 1:10
        psi(j) = exp(-1*(x(i)-c(j))^2 /(2*var(j))); end
    for j = 1:10
        phi(j) = psi(j)*x(i)/sum(psi); end
    f = phi * w;
    plotvar_f = [plotvar_f;psi];
    zdot(i+1) = alpha_z*(beta_z*(g-y(i))-ydot(i))+f;
    ydot(i+1) = ydot(i) + zdot(i)*dt;
    y(i+1) = y(i) + ydot(i+1)*dt;
    x(i+1) = x(i) - alpha_x*x(i)*dt;
end
 end





