n = 79;
f = zeros(n^2, 1);
% Boundary pts.
b = zeros(n^2, 1);
% Width of evenly spaced mesh.
dx = 1 / (n+1);
% Height of evenly spaced mesh.
dy = 1/ (n+1);
for j = 1:n
    for i = 1:n
        f((j-1) * n + i) = 2 / (1 + i * dx)^3 + 2 / (1 + j * dy)^3;
    end
end
% u evaluate at each internal grid pt.
u_tilde = zeros(n^2, 1);
for j = 1:n
    for i = 1:n
        u_tilde((j-1) * n + i) = 1 / (1 + i * dx) + 1 / (1 + j * dy);
    end
end
% The coefficients of u_ij loaded in by row.
A = sparse(n^2, n^2);
I = zeros(n, n);

for j = 1:n
    for i = 1:n
        I(i,j) = (j-1) * n + i;
    end
end

C = [1 / dx^2, -2 / dx^2 - 2 / dy^2, 1 / dx^2, 1 / dy^2, 1 / dy^2];

for j = 1:n
    for i = 1:n
        temp = 0;
        A((j-1) * n + i, I(i,j)) = C(2);
        if i ~= n
            A((j-1) * n + i, I(i+1,j)) = C(1);
        else
            A((j-1) * n + i, I(i,j)) = A((j-1) * n + i, I(i,j)) + C(1);
            temp = temp - C(1) * dx / 4;
        end
        if i ~= 1
            A((j-1) * n + i, I(i-1,j)) = C(3);
        else
            temp = temp + C(3) * (1 + 1 / (1 + j * dy));
        end
        if j ~= n
            A((j-1) * n + i, I(i, j+1)) = C(4);
        else
            temp = temp + C(4) * (1 / 2 + 1 / (1 + i * dx));
        end
        if j ~= 1
            A((j-1) * n + i, I(i, j-1)) = C(5);
        else
            temp = temp + C(5) * (1 + 1 / (1 + i * dx));
        end
        
        b((j-1) * n + i) = temp;
    end
end
u = A \ (f - b);
result1 = max(abs(u - u_tilde));
u_bdy = zeros(n, 1);
u_bdy_tilde = zeros(n, 1);
for j = 1:n
    u_bdy(j) = u(j * n) - dx / 4;
    u_bdy_tilde(j) = 1 / (1 + (n+1) * dx) + 1 / (1 + j * dy);
end
result2 = max(abs(u_bdy - u_bdy_tilde));
result = max([result1;result2]);
s = 'The maximum error in the numerical solution for n = %d is: %1.8f.';
fprintf(s, n, result);