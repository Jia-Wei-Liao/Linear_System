function [x, ii, rho, re] = GaussSeidel(A, b)
D = diag(diag(A));
R = -triu(A, 1);
L = -tril(A, -1);

M = D-L;
N = R;
H = M\N;
f = M\b;
rho = max(abs(eig(H)));

n = size(b, 1);
x = rand(n, 1);
ii = 0;
tol = 1e-16;
re = 1;
MaxStep = 10000000;

if rho >= 1
    x = NaN;
    re = NaN;
    fprintf('The iteration of Gauss Seidel does not converge.\n');
    
else
    while re > tol & ii < MaxStep
        ii = ii + 1;
        x_ = H*x + f;
        re = RelativeError(x_, x);
        x = x_;
    end
    fprintf('The number of Gauss Seidel iteration is %d.\n', ii);
end

end