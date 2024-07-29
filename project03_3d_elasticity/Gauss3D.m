function [xi, eta, zeta, w] = Gauss3D(N1, N2, N3)

% preallocation
xi = zeros(N1*N2*N3,1);
eta = xi;
zeta=xi;
w = xi;

% generate 1D rule
[x1, w1] = Gauss(N1, -1, 1);

[x2, w2] = Gauss(N2, -1, 1);

[x3, w3] = Gauss(N3, -1, 1);

for ii = 1 : N1
    for jj = 1 : N2
        for kk = 1 : N3
            xi(  (ii-1)*N1*N2 + (jj-1)*N1 + kk ) = x1(ii);
            eta( (ii-1)*N1*N2 + (jj-1)*N1 + kk ) = x2(jj);
            zeta((ii-1)*N1*N2 + (jj-1)*N1 + kk ) = x3(kk);
            w(   (ii-1)*N1*N2 + (jj-1)*N1 + kk ) = w1(ii) * w2(jj) * w3(kk);
        end
    end
end

% EOF