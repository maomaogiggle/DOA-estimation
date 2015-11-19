function x0 = linear_search(W,dxx,dx,K,N,delta)
id = 1:N;
V = nchoosek(id,K);
num = nchoosek(N,K);
e = ones(K,1);
%searching
y = zeros(num,1);
for i = 1:num
    if (e'*dx(V(i,:))==0)
        flag = 0;
        for k = 1:size(W,2)
            v = W(:,k);
            if ((1/K)*abs(v(V(i,:))'*e) <= delta)
                flag = flag + 1;
            end
        end
        if (flag == size(W,2))
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



