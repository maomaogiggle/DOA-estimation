function x0 = linear_amb(V,dxx,dx,K,N,delta)
e = ones(N,1);
mu = 10;
Q = 4;
x0 = [ones(K/2,1);zeros(N-K,1);ones(K/2,1)];
for q = 1:Q
    cvx_begin
        variable x(N);
        minimize (-x'*dxx+mu*(e'*x-2*x0'*x));
        subject to
        0 <= x'*dx <= 0;
        K <= e'*x <= K;
        0 <= x <= 1;
        for i = 1:size(V,2)
            W = (1/(K*K))*real(V(:,i)*V(:,i)');
            x'*W*x <= delta;
        end
    cvx_end
    x0 = x;
end
end