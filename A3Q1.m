m = 10;
% Boundary pts.
b = zeros(m, 1);
% Length of evenly spaced mesh.
h = 1 / m;

c = 1/2 * (exp(1) - exp(-1)) - pi;

% y evaluate at each internal grid pt.
y = zeros(m, 1);

for j = 1:m
    x = j * h;
    y(j) = 1 / 2 * (exp(x) + exp(-x)) + sin(pi * x);
end
% The coefficients of a_ij loaded in by row.
A = sparse(m, m);

for i = 1:m
    for j = 1:m
        if i - j == 0
            A(i,j) = 2 / h + 2 * h / 3;
        elseif i - j == 1
            A(i,j) = -1 / h + h * (1 / 6);
        elseif i - j == -1
            A(i,j) = -1 / h + h * (1 / 6);
        end
    end
    if i ~= m
        b(i) = -h/6*((i-1)*h*c + 4*i*h*c + 6 + (i+1)*h*c) + ...
        ((1 + pi^2) / (pi^2*h))*(2*sin(pi*i*h) ...
        - sin(pi*(i-1)*h) - sin(pi*(i+1)*h));
    else
        b(i) = -h/6*((i-1)*h*c + 2*i*h*c + 3) + ...
        ((1 + pi^2) / (pi^2*h))*(sin(pi*i*h) ...
        - sin(pi*(i-1)*h) - pi*h*cos(pi*i*h));
    end

end
y_m_tilde = A \ b;
y_m = zeros(m, 1);

for i = 1:m
    x = i * h;
    y_m(i) = 1 + c*x + y_m_tilde(i);
end
result = max(abs(y_m - y));
s = 'The maximum error in the numerical solution for m = %d is: %1.8f.';
fprintf(s, m, result);