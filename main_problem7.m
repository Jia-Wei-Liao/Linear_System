clc; clear; close all;
addpath(genpath('function'));

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
R = -triu(A, 1);
L = -tril(A, -1);

% SOR for best w
wb = 2/(1+sqrt(1-max(abs(eig(L+R)))^2));
[x, best_ii, rho, re] = SOR(A, b, wb);

% SOR for general 1<w<2
N = 2e+4;
k = 1;
dw = 1/N;
w0 = 0+dw; w_end = 2-dw;
I = eye(4);
Record = zeros(N, 3);

for w = w0:dw:w_end
    [x, iter, ~, re] = SOR(A, b, w);
    Lw = (I-w*L)\((1-w)*I+w*R);
    rho = max(abs(eig(Lw)));
    Record(k, :) = [w iter, rho];
    k = k + 1;
end
Record = Record(Record(:,1)>0, :);

% Plot the Figure 1.
figure(1)
semilogy(Record(:, 1), Record(:, 2), 'LineWidth', 2);
xlabel('$w$', 'fontsize', 12, 'interpreter', 'latex');
ylabel('iteration $i$', 'fontSize', 12, 'interpreter', 'latex');
grid on
hold on
semilogy(wb, best_ii, '-o', 'MarkerFaceColor', 'r');
hold off

% Plot the Figure 2.
figure(2)
plot(Record(:, 1), Record(:, 3), '.-');
xlabel('$w$', 'interpreter', 'latex', 'fontsize', 12);
ylabel('$\rho(L_\omega)$', 'interpreter', 'latex', 'fontsize', 12);
ylim([-0.3, 1]);
grid on
hold on
[~, index] = min(Record(:, 3));
plot(Record(index, 1), Record(index, 3), 'ro', 'MarkerFaceColor', 'red');
