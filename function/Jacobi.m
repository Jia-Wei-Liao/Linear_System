function [x, ii, rho, re] = Jacobi(A, b)
D = diag(diag(A));
R = -triu(A, 1);
L = -tril(A, -1);

M = D;
N = L+R;
J = M\N;
f = M\b;
rho = max(abs(eig(J)));

n = size(b, 1);
x = rand(n, 1);
ii = 0;
tol = 1e-16;
re = 1;
MaxStep = 10000000;

if rho >= 1
    x = NaN;
    re = NaN;
    fprintf('The iteration of Jacobi method does not converge.\n');
    
else
    while re > tol & ii < MaxStep
        ii = ii + 1;
        x_ = J*x + f;
        re = RelativeError(x_, x);
        x = x_;
    end
    fprintf('The number of Jacobi iteration is %d.\n', ii);
end

end