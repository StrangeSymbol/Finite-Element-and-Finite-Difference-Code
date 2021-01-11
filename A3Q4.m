m = 9;
% Boundary pts.
b = zeros(m^2, 1);
% Width of evenly spaced mesh.
h = 1 / (m+1);

for l = 1:m
    for k = 1:m
        c1 = -8*h/3*(4*(k*h)^2 - 8*k*h + (k*h + (k-1)*h)^2 - 2*(k-1)*h ...
             - 2*(k+1)*h + (k*h + (k+1)*h)^2);
        c2 = -8*h/3*(4*(l*h)^2 - 8*l*h + (l*h + (l-1)*h)^2 - 2*(l-1)*h ...
             - 2*(l+1)*h + (l*h + (l+1)*h)^2);
        b((l-1) * m + k) = c1*h + c2*h;
    end
end
% y evaluate at each internal grid pt.
u = zeros(m^2, 1);

for j = 1:m
    for i = 1:m
        x = i*h;
        y = j*h;
        u((j-1) * m + i) = 16*x*(1-x)*y*(1-y);
    end
end
% The coefficients of a_ij loaded in by row.
A = sparse(m^2, m^2);

for j = 1:m
    for i = 1:m
        for l = 1:m
            for k = 1:m
                if abs(j-l) == 1 && abs(i-k) <= 1
                   A((j-1)*m+i,(l-1)*m+k) = -1 / 3;
                elseif abs(j-l) == 0 && abs(i-k) == 0
                   A((j-1)*m+i,(l-1)*m+k) = 8/3;
                elseif abs(j-l) <= 1 && abs(i-k) == 1
                   A((j-1)*m+i,(l-1)*m+k) = -1 / 3;
                end
            end
        end
    end
end
u_m = A \ b;
result = max(abs(u_m - u));
s = 'The maximum error in the numerical solution for m = %d is: %1.8f.';
fprintf(s, m, result);