clear all; % clean memory
clc; % clean the screen

n_int = 3;

[xi, weight] = Gauss(n_int, -1, 1);

n = 6;

exact = (1 - (-1)^(n+1)) / (n+1);

appro = 0;
for l = 1 : n_int
    appro = appro + weight(l) * xi(l)^n;
end

error = exact - appro