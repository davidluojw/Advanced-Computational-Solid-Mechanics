% This is the Newton-Raphson algorithm solver for HW3
clear all; clc;

% exact internal force
x = 15;
N1 = @(d1,d2) x * d1 / (10.0 - d1) - 0.5 * d2*d2;
N2 = @(d1,d2) d2 - d1;
load_step = 0.25;
load_max = 10.0;
load_num = load_max / load_step;


% external force
F1 = @(s) 0.25 * s;
F2 = @(s) 0.0;

% Jacobian
J = zeros(2,2);
dN1_dd1 = @(d1,d2) 10*x / ((10.0-d1)*(10.0-d1));
dN1_dd2 = @(d1,d2) -d2;
dN2_dd1 = @(d1,d2) -1.0;
dN2_dd2 = @(d1,d2) 1.0;

% Newton-Raphson Algorithm

for i = 0:load_num
    

