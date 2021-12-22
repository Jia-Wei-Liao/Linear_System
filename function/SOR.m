function [x, ii, rho, re] = SOR(A, b, wb)
D = diag(diag(A));
R = -triu(A, 1);
L = -tril(A, -1);

Mw = D-wb*L;
Nw = (1-wb)*D + wb*R;
S = Mw\Nw;
f = (Mw\b)*wb;
rho = max(abs(eig(S)));

n = size(b, 1);
x = zeros(n, 1); %rand(n, 1);
ii = 0;
re = 1;
tol = 1e-16;
MaxStep = 10000000;

if rho >= 1
    x = NaN;
    re = NaN;
    fprintf('The iteration of SOR does not converge.\n');
    
else
    while re > tol & ii < MaxStep
        ii = ii + 1;
        x_ = S*x + f;
        re = RelativeError(x_, x);
        x = x_;
    end
    fprintf('The number of SOR iteration is %d.\n', ii);
end

end