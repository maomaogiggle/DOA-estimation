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
P_o = zeros(K,2,length(phi_s));


for i = 1:length(phi_s)
    
    alpha = (sin(phi_s(i))^2)/(cos(theta_s)^2)+(cos(phi_s(i))^2)/(sin(theta_s)^2);
    beta = (cos(phi_s(i))^2)/(cos(theta_s)^2)+(sin(phi_s(i))^2)/(sin(theta_s)^2);
    zeta = (sin(2*phi_s(i)))/(sin(theta_s)^2)-(sin(2*phi_s(i)))/(cos(theta_s)^2);
    x_o = directive_search(dxx,dx,dyy,dy,dxy,K,N1*N2,alpha,beta,zeta);

    index = 1;
    for k = 1:N1*N2
        if (abs(x_o(k)-1)<=0.1)
            P_o(index,:,i) = mn(k,:);
            index = index + 1;
        end
    end
end

save('P_o.mat','P_o');








