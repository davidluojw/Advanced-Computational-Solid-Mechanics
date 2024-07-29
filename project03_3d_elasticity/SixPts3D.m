function [xi, eta, zeta, w] = SixPts3D()

% preallocation
xi = zeros(6,1);
eta = xi;
zeta = xi;
w = 4 / 3 * ones(6,1);


xi(1) = 1; xi(2) = -1;
eta(3) = 1; eta(4) = -1;
zeta(5) = 1; zeta(6) = -1;


% EOF