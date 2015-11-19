function x0 = planar_search_n(dxx,dx,dyy,dy,dxy,K,N)
id = 1:N;
V = nchoosek(id,K);
num = nchoosek(N,K);
e = ones(K,1);
%searching
y = zeros(num,1);
for i = 1:num
    if ((e'*dx(V(i,:))==0) && (e'*dy(V(i,:))==0))
        if ((e'*dxy(V(i,:))==0) && (e'*dxx(V(i,:)) == e'*dyy(V(i,:))))
           y(i) = e'*dxx(V(i,:));
        end
   end  
end
[ymax,Imax] = max(y);
x0 = zeros(N,1);
for i = 1:K
x0(V(Imax,i)) = 1;
end
end



