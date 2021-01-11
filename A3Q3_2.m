m = 9;
lambda = 100;
% Boundary pts.
b = zeros(m, 1);
% The number of points on far left and right each.
num_rapid = floor(m / 6);
% Number of points in the middle.
num_middle = ceil(2*m / 3) + 1;
if mod(ceil(2*m / 3), 2) ~= 0
    num_middle = ceil(2*m / 3);
end
% Length of spaced mesh.
h = 0.92 / (num_middle+1);
p = 0.04 / num_rapid; % The equal spacing for left and right side.
Q_0 = 1;

x = zeros(m+2,1);

for i = 2:m+2
    if i <= num_rapid + 1
        x(i) = (i-1) * p;
    elseif i <= num_rapid + num_middle + 2
        x(i) = h + x(i-1);
    else
        x(i) = p + x(i-1);
    end
end

% y evaluate at each internal grid pt.
y = zeros(m, 1);
c2 = (exp(lambda) - 1) / (exp(lambda) - exp(-lambda));
c1 = (1 - exp(-lambda)) / (exp(lambda) - exp(-lambda));
for j = 1:m
    y(j) = c1 * exp(lambda*x(j+1)) + c2 * exp(-lambda*x(j+1));
end
% The coefficients of a_ij loaded in by row.
A = sparse(m, m);

for i = 1:m
    for j = 1:m
        d1 = x(i+1) - x(i);
        d2 = x(i+2) - x(i+1);
        if i - j == 0
            A(i,j) = 1/d1 + 1/d2 + lambda^2*(d1 / 3 + d2 / 3);
        elseif i - j == 1
            A(i,j) = -1 / d1 + lambda^2*d1 / 6;
        elseif i - j == -1
            A(i,j) = -1 / d2 + lambda^2*d2 / 6;
        end
    end
    b(i) = 1/2*(-lambda^2 * d1 - lambda^2 * d2);
end
y_m = Q_0 + A \ b;
result = max(abs(y_m - y));
s = 'The maximum error in the numerical solution for m = %d is: %1.8f.';
fprintf(s, m, result);