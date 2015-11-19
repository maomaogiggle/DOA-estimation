function w0 = dinkelbach(dxx,dx,dyy,dy,dxy,K,N,alpha,beta,zeta,eta)

%compose matrix A and B
A = (alpha/K)*dxx*ones(1,N)+(beta/K)*dyy*ones(1,N) + (zeta/K)*dxy*ones(1,N);
B = dxx*dyy'-dxy*dxy';
F = 1;
kesi = 0.01;
ee = ones(N,1);

while(abs(F) > kesi)
    cvx_begin
        variable W(N,N); variable w(N,1);
        minimize (trace(W*(A-eta*B)));
        subject to
            0 <= w'*dx <= 0;
            0 <= w'*dy <= 0;
            K <= w'*ee <= K;
            [W,w;w',1] == semidefinite(N+1);
            for k = 1:N
                e = [zeros(k-1,1);ones(1,1);zeros(N-k,1)];
                E = e*e';
                0 <= trace(W*E)-e'*w <= 0;
            end
    cvx_end
    F = trace(W*(A-eta*B));
    eta = trace(W*A)/trace(W*B);
    
end
[Ys,Is] = sort(w,'descend');
w0 = zeros(N,1);
w0(Is(1:K))=1;
end
