clc; clear; close all;

A = [
        3 -0.25 -0.25    0 ;
    -0.25     5     0 -0.5 ;
    -0.25     0     2    0 ;
        0 -0.25     0    3 ;
];
b = rand(4, 1);
D = diag(diag(A));

A = D\A;
b = D\b;

% SOR
R = -triu(A, 1);
L = -tril(A, -1);
wb = 2/(1+sqrt(1-max(abs(eig(L+R)))^2));
[x, ii, rho, re] = SOR(A, b, wb);

for w = 0.1:0.1:1.9
    [x, ii, rho, re] = SOR(A, b, w);
end


