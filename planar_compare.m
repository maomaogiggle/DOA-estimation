%optimum subarray versus the elevation angle
clear;
clc;

N1 = 5;
N2 = 5;
K = 10;
lambda = 1;
d = lambda/2;
k0 = 2*pi/lambda;
n = -(N1-1)/2:1:(N1-1)/2;
m = [-(N2-1)/2:1:(N2-1)/2]';
mn = [];
for i = 1:N1
    mn = [mn;[n(i)*ones(N2,1),m]]; %mn is the position vector of N*N/2 by 2 dimension
end

%define the array position matrix
dx = d*mn(:,1);
dy = d*mn(:,2);
dxx = dx.*dx;
dyy = dy.*dy;
dxy = dx.*dy;

%generate the arraival direction
theta_s = 10*pi/180;
phi_s = [0:1:180]*pi/180;

snr = 10; %signal to noise ratio

%read the optimum configuration corresponding to each azimuth angle
load('P_1.mat');

e = ones(1,K);
rho = 10^(snr/10);

mse_theta_ot = zeros(length(phi_s),1);
mse_phi_ot = zeros(length(phi_s),1);

for i = 1:length(phi_s)
    
    P_o = P_1(:,:,i);

    mse_theta_ot(i) = (1+K*rho)./(2*K*(rho.^2)*(4*pi*pi/lambda/lambda)*(cos(theta_s).^2))...
             *(d*d*e*(P_o(:,1).^2)*(sin(phi_s(i)).^2)+d*d*e*(P_o(:,2).^2)*(cos(phi_s(i)).^2)...
             -d*d*e*(P_o(:,1).*P_o(:,2))*sin(2*phi_s(i)))/(d*d*e*(P_o(:,1).^2)*d*d*e*(P_o(:,2).^2)...
             -d*d*e*(P_o(:,1).*P_o(:,2))*d*d*e*(P_o(:,1).*P_o(:,2)));
    mse_phi_ot(i) = (1+K*rho)./(2*K*(rho.^2)*(4*pi*pi/lambda/lambda)*(sin(theta_s).^2))...
             *(d*d*e*(P_o(:,1).^2)*(cos(phi_s(i)).^2)+d*d*e*(P_o(:,2).^2)*(sin(phi_s(i)).^2)...
             +d*d*e*(P_o(:,1).*P_o(:,2))*sin(2*phi_s(i)))/(d*d*e*(P_o(:,1).^2)*d*d*e*(P_o(:,2).^2)...
             -d*d*e*(P_o(:,1).*P_o(:,2))*d*d*e*(P_o(:,1).*P_o(:,2)));
end
plot(phi_s*180/pi,10*log10(mse_phi_ot+mse_theta_ot),'c');










