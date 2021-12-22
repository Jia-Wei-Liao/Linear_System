clc; clear; close all;
addpath(genpath('function'));

n = 5;     % 5, 10
wb = 1.1;  % 1.1, 1.2, 1.8

A = hilb(n);
b = ones(n, 1);

%% Jacobi
[x, ii, rho, re] = Jacobi(A, b);

%% Gauss Seidel
[x, ii, rho, re] = GaussSeidel(A, b);

%% SOR
[x, ii, rho, re] = SOR(A, b, wb);

