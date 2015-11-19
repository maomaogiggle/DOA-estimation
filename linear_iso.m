%optimum isotorpic array
clear;
clc;

N = 20;
K = 10;
n = -(N-1)/2:1:(N-1)/2;
id = 1:N;
lambda = 1;
d = lambda/2;
k0 = 2*pi/lambda;

%define the array position matrix
dx = d*n';
dxx = dx.*dx;

%avoiding ambiguity
vjs1 = exp(1i*k0*d*n'*0.12);
vjs2 = exp(1i*k0*d*n'*0.24);
vjs3 = exp(1i*k0*d*n'*0.36);
vjs4 = exp(1i*k0*d*n'*0.48);
W = [vjs1,vjs2,vjs3,vjs4];
delta = 0.2;
% x_am = linear_amb(W,dxx,dx,K,N,delta);
x_am = linear_iso_search(W,dxx,dx,K,N,delta);
%calculate the position
P_am = zeros(K,1);
index = 1;
 for i = 1:N
     if (abs(x_am(i)-1)<=0.4)
         P_am(index) = n(i);
         index = index + 1;
     end
 end
%optimum isotropic array without considering ambiguity
x_o = [ones(K/2,1);zeros(N-K,1);ones(K/2,1)];
P_o = zeros(K,1);
index = 1;
 for i = 1:N
     if (abs(x_o(i)-1)<=0.4)
         P_o(index) = n(i);
         index = index + 1;
     end
 end

%compare beampatterns
vs_o = exp(-1i*k0*d*P_o*0);
vs_am = exp(-1i*k0*d*P_am*0);
theta = 0:0.5:180;
u = cos(theta*pi/180);
v = exp(1i*k0*d*P_o*u);
figure;
plot(theta,20*log10(abs(vs_o'*v)/max(abs(vs_o'*v))));
hold on;
v = exp(1i*k0*d*P_am*u);
plot(theta,20*log10(abs(vs_am'*v)/max(abs(vs_am'*v))),'k');

%calculate the estimation variance
snr = -5:1:30; %signal to noise ratio
T = 1; %the number of snapshots
theta_search = 0:0.0005:pi;
V_o = exp(1i*k0*d*P_o*cos(theta_search));
V_am = exp(1i*k0*d*P_am*cos(theta_search));
V_ula = exp(1i*k0*d*n'*cos(theta_search));
Q = 500;%Q Monte-Carlo runs

mse_o = zeros(length(snr),1);
mse_am = zeros(length(snr),1);
mse_ula = zeros(length(snr),1);

%generate the arraival direction
theta_s  =  pi/3;
vs_o = exp(1i*k0*d*P_o*cos(theta_s));
vs_am = exp(1i*k0*d*P_am*cos(theta_s));
%ULA
vs_ula = exp(1i*k0*d*n'*cos(theta_s));

for q = 1:length(snr)
    variance_o = 0;
    variance_am = 0;
    variance_ula = 0;
    for i = 1:Q
        %generate the estimated signal
        X = (randn(1,T) + 1i * randn(1,T)) / sqrt(2);
        Y_o = vs_o * X;
        Y_am = vs_am * X;
        Y_ula = vs_ula * X;
        % generate the white noise
        sigma_o = 10^(-snr(q)/10) * (norm(Y_o,'fro')^2) / (K * T);
        sigma_am = 10^(-snr(q)/10) * (norm(Y_am,'fro')^2) / (K * T);      
        sigma_ula = 10^(-snr(q)/10) * (norm(Y_ula,'fro')^2) / (N * T);           
        % noise
        E_o = sqrt(sigma_o)/sqrt(2)*randn(K,T) + 1i*sqrt(sigma_o)/sqrt(2)*randn(K,T);
        E_am = sqrt(sigma_am)/sqrt(2)*randn(K,T) + 1i*sqrt(sigma_am)/sqrt(2)*randn(K,T);
        E_ula = sqrt(sigma_ula)/sqrt(2)*randn(N,T) + 1i*sqrt(sigma_ula)/sqrt(2)*randn(N,T);        
        %received signal
        Y_o = Y_o + E_o;
        Y_am = Y_am + E_am;   
        Y_ula = Y_ula + E_ula;           
        %DOA estimation
        y_o = abs(V_o'*Y_o);
        y_am = abs(V_am'*Y_am);  
        y_ula = abs(V_ula'*Y_ula);          
        [y_m1,I_max] = max(y_o);
        theta_est_o = theta_search(I_max);
        [y_m2,I_max] = max(y_am);
        theta_est_am = theta_search(I_max);
        [y_m3,I_max] = max(y_ula);
        theta_est_ula = theta_search(I_max);
        %variance
        variance_o = variance_o + (theta_est_o-theta_s)^2;
        variance_am = variance_am + (theta_est_am-theta_s)^2;    
        variance_ula = variance_ula + (theta_est_ula-theta_s)^2;            
    end
    mse_o(q) = variance_o/Q;
    mse_am(q) = variance_am/Q; 
    mse_ula(q) = variance_ula/Q; 
end

%array without ambiguity solved
e = ones(1,K);
figure;
plot(snr,10*log10(mse_o));
hold on;
rho = 10.^(snr/10);
mse_ot= (1+K*rho)./(2*K*(rho.^2)*(4*pi*pi/lambda/lambda)*(sin(theta_s)^2)*d*d*(e*(P_o.^2)));
plot(snr,10*log10(mse_ot),'k');

%array with ambiguity solved
hold on;
plot(snr,10*log10(mse_am),'g');
hold on;
mse_amt= (1+K*rho)./(2*K*(rho.^2)*(4*pi*pi/lambda/lambda)*(sin(theta_s)^2)*d*d*(e*(P_am.^2)));
plot(snr,10*log10(mse_amt),'y');

%ULA
hold on;
plot(snr,10*log10(mse_ula),'r');
hold on;
mse_ulat= (1+N*rho)./(2*N*(rho.^2)*(4*pi*pi/lambda/lambda)*(sin(theta_s)^2)*d*d*((n.^2)*ones(N,1)));
plot(snr,10*log10(mse_ulat),'c');

