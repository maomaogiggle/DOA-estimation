function x0 = directive_search(dxx,dx,dyy,dy,dxy,K,N,alpha,beta,zeta)
id = 1:N;
V = nchoosek(id,K);
num = nchoosek(N,K);
e = ones(K,1);
%searching
y = zeros(num,1);
for i = 1:num
    if ((e'*dx(V(i,:))==0) && (e'*dy(V(i,:))==0))
        y(i) = (e'*dxx(V(i,:))*e'*dyy(V(i,:))-e'*dxy(V(i,:))*e'*dxy(V(i,:)))...
            /(alpha*e'*dxx(V(i,:))+beta*e'*dyy(V(i,:))+zeta*e'*dxy(V(i,:)));
   end  
end
[ymax,Imax] = max(y);
x0 = zeros(N,1);
for i = 1:K
x0(V(Imax,i)) = 1;
end
end

