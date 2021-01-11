m = 19;
lambda = sqrt(10);
% Boundary pts.
b = zeros(m, 1);
% Length of evenly spaced mesh.
h = 1 / (m+1);
Q_0 = 1;
% y evaluate at each internal grid pt.
y = zeros(m, 1);
c2 = (exp(lambda) + 1 / lambda^2 - 1) / (exp(lambda) - exp(-lambda));
c1 = 1 - c2;
for j = 1:m
    x = j * h;
    y(j) = c1 * exp(lambda*x) + c2 * exp(-lambda*x) + x / lambda^2;
end
% The coefficients of a_ij loaded in by row.
A = sparse(m, m);

for i = 1:m
    for j = 1:m
        if i - j == 0
            A(i,j) = 2 / h + 2 * h / 3 * lambda ^ 2;
        elseif i - j == 1
            A(i,j) = -1 / h + lambda^2 * h * (1 / 6);
        elseif i - j == -1
            A(i,j) = -1 / h + lambda^2 * h * (1 / 6);
        end
    end
    b(i) = i * h^2 - lambda^2 * h;
end
y_m = 1 + A \ b;
result = max(abs(y_m - y));
s = 'The maximum error in the numerical solution for m = %d is: %1.8f.';
fprintf(s, m, result);