%optimum isotorpic planar array
clear;
clc;

N1 = 5;
N2 = 5;
K = 10;
n = -(N1-1)/2:1:(N1-1)/2;
m = [-(N2-1)/2:1:(N2-1)/2]';
mn = [];
for i = 1:N1
    mn = [mn;[n(i)*ones(N2,1),m]]; %mn is the position vector of N*N/2 by 2 dimension
end
lambda = 1;
d = lambda/2;
k0 = 2*pi/lambda;

%define the array position matrix
dx = d*mn(:,1);
dy = d*mn(:,2);
dxx = dx.*dx;
dyy = dy.*dy;
dxy = dx.*dy;

%optimum isotropic array with ambiguity
x_o = planar_search_n(dxx,dx,dyy,dy,dxy,K,N1*N2);
P_o = zeros(K,2);
index = 1;
 for i = 1:N1*N2
     if (abs(x_o(i)-1)<=0.4)
         P_o(index,:) = mn(i,:);
         index = index + 1;
     end
 end
% %plot the selected subarray
% spy(fliplr(reshape(x_o,5,5)'));

%generate the arraival direction
theta_s = 10*pi/180;
phi_s = 175*pi/180;
ux_s  =  sin(theta_s)*cos(phi_s);
uy_s = sin(theta_s)*sin(phi_s);
vs_o = exp(1i*k0*d*P_o*[ux_s;uy_s]);
%compare beampatterns
theta = 0:0.02:pi/2;
phi = 0:0.02:2*pi;
Y_o = zeros(length(theta),length(phi));
for i = 1:length(theta)
    for k = 1:length(phi)
        ux = sin(theta(i))*cos(phi(k));
        uy = sin(theta(i))*sin(phi(k));
        v = exp(1i*k0*d*P_o*[ux;uy]);
        Y_o(i,k) = abs(vs_o'*v);
    end
end
[X,Y] = meshgrid(phi,theta);
mesh(X,Y,Y_o);

%calculate the estimation variance
snr = -5:1:30; %signal to noise ratio
T = 1; %the number of snapshots
Q = 500;%Q Monte-Carlo runs
%searching grid points
theta = 0:0.002:pi/2;
phi = 0:0.002:2*pi;
ux = sin(theta')*cos(phi);
uy = sin(theta')*sin(phi);
u_search = [ux(:)';uy(:)'];

mse_theta_o = zeros(length(snr),1);
mse_phi_o = zeros(length(snr),1);

for q = 1:length(snr)
    variance_theta_o = 0;
    variance_phi_o = 0;
    for i = 1:Q
        %generate the estimated signal
        X = (randn(1,T) + 1i * randn(1,T)) / sqrt(2);
        Y_o = vs_o * X;
        % generate the white noise
        sigma_o = 10^(-snr(q)/10) * (norm(Y_o,'fro')^2) / (K * T);      
        % noise
        E_o = sqrt(sigma_o)/sqrt(2)*randn(K,T) + 1i*sqrt(sigma_o)/sqrt(2)*randn(K,T);
        %received signal
        Y_o = Y_o + E_o;
        %DOA estimation
        v_o = exp(1i*k0*d*P_o*u_search);
        y_o = abs(v_o'*Y_o);
        %searching the peak
        [y_m1,I_max] = max(y_o);
        if (mod(I_max,length(theta)) == 0)
            theta_est_o = theta(length(theta));
        else
            theta_est_o = theta(mod(I_max,length(theta)));
        end
        phi_est_o = phi(floor(I_max/length(theta))+1);
        %variance
        variance_theta_o = variance_theta_o + (theta_est_o-theta_s)^2;
        variance_phi_o = variance_phi_o + (phi_est_o-phi_s)^2;      
    end
    mse_theta_o(q) = variance_theta_o/Q;
    mse_phi_o(q) = variance_phi_o/Q;
end

%array without ambiguity solved
e = ones(1,K);
figure;
plot(snr,10*log10(mse_theta_o));
hold on;
rho = 10.^(snr/10);
mse_theta_ot= (1+K*rho)./(2*K*(rho.^2)*(4*pi*pi/lambda/lambda)*(cos(theta_s)^2)...
             *((cos(phi_s)^2)*d*d*e*(P_o(:,1).^2)+(sin(phi_s)^2)*d*d*e*(P_o(:,2).^2)));
plot(snr,10*log10(mse_theta_ot),'k');
hold on;
plot(snr,10*log10(mse_phi_o),'r');
hold on;
mse_phi_ot= (1+K*rho)./(2*K*(rho.^2)*(4*pi*pi/lambda/lambda)*(sin(theta_s)^2)...
             *((sin(phi_s)^2)*d*d*e*(P_o(:,1).^2)+(cos(phi_s)^2)*d*d*e*(P_o(:,2).^2)));
plot(snr,10*log10(mse_phi_ot),'c');




