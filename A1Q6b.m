n = 79;
% Width of evenly spaced mesh.
dx = 1 / (n+1);
% Height of evenly spaced mesh.
dy = 1 / (n+1);

count = 0;
P = [];
for j = -n:n
    for i = -n:n
        if (i * dx)^2 + (j * dy)^2 < 1
            P = [P;[i,j]];
            count = count + 1;
        end
    end
end
f = zeros(count, 1);
% Boundary pts.
b = zeros(count, 1);

for i = 1:count
    f(i) = 16 * ((P(i,1) * dx)^2 + (P(i,2) * dy)^2);
end
% u evaluate at each internal grid pt.
u_tilde = zeros(count, 1);

for i = 1:count
    u_tilde(i) = ((P(i,1) * dx)^2 + (P(i,2) * dy)^2)^2;
end
% The coefficients of u_ij loaded in by row.
A = sparse(count, count);

C = [1 / dx^2, -2 / dx^2 - 2 / dy^2, 1 / dx^2, 1 / dy^2, 1 / dy^2];

for m = 1:count
    i = P(m,1);
    j = P(m,2);
    temp = 0;
    % Interior pts.
    if ((i+1) * dx)^2 + (j * dy)^2 < 1 && ((i-1) * dx)^2 + (j * dy)^2 < 1 
        A(m, getIndex(P,i+1,j,count)) = C(1);
        A(m, getIndex(P,i,j,count)) = -2 / dx^2;
        A(m, getIndex(P,i-1,j,count)) = C(3);
    else
        if ((i+1) * dx)^2 + (j * dy)^2 > 1
            tau = abs(sqrt(1 - (j*dy)^2) - i*dx) / dx;
            A(m, getIndex(P,i,j,count)) =  -(3 - tau) / (tau * dx^2);
            A(m, getIndex(P,i-1,j,count)) = 2 * (2 - tau) / ((tau + 1) * dx^2);
            A(m, getIndex(P,i-2,j,count)) = (tau - 1) / ((tau + 2) * dx^2);
            temp = temp + 6 / ((tau * (tau + 1) * (tau + 2)) * dx^2);
        elseif ((i+1) * dx)^2 + (j * dy)^2 == 1.0
            A(m, getIndex(P,i,j,count)) = -2 / dx^2;
            A(m, getIndex(P,i-1,j,count)) = C(3);
            temp = temp + C(1);
        end
        if ((i-1) * dx)^2 + (j * dy)^2 > 1
            tau = abs(-sqrt(1 - (j*dy)^2) - i*dx) / dx;
            A(m, getIndex(P,i,j,count)) =  -(3 - tau) / (tau * dx^2);
            A(m, getIndex(P,i+1,j,count)) = 2 * (2 - tau) / ((tau + 1) * dx^2);
            A(m, getIndex(P,i+2,j,count)) = (tau - 1) / ((tau + 2) * dx^2);
            temp = temp + 6 / ((tau * (tau + 1) * (tau + 2)) * dx^2);
        elseif ((i-1) * dx)^2 + (j * dy)^2 == 1.0
            A(m, getIndex(P,i+1,j,count)) = C(1);
            A(m, getIndex(P,i,j,count)) = -2 / dx^2;
            temp = temp + C(3);
        end
    end
    if (i * dx)^2 + ((j+1) * dy)^2 < 1 && (i * dx)^2 + ((j-1) * dy)^2 < 1
        A(m, getIndex(P,i,j,count)) = A(m, getIndex(P,i,j,count)) - 2 / dy^2;
        A(m, getIndex(P,i,j+1,count)) = C(4);
        A(m, getIndex(P,i,j-1,count)) = C(5);
    else
        if (i * dx)^2 + ((j+1) * dy)^2 > 1
            tau = abs(sqrt(1 - (i*dx)^2) - j*dy) / dy;
            A(m, getIndex(P,i,j,count)) = A(m, getIndex(P,i,j,count)) - (3 - tau) / (tau * dy^2);
            A(m, getIndex(P,i,j-1,count)) = 2 * (2 - tau) / ((tau + 1) * dy^2);
            A(m, getIndex(P,i,j-2,count)) = (tau - 1) / ((tau + 2) * dy^2);
            temp = temp + 6 / ((tau * (tau + 1) * (tau + 2)) * dy^2);
        elseif (i * dx)^2 + ((j+1) * dy)^2 == 1.0
            A(m, getIndex(P,i,j,count)) = A(m, getIndex(P,i,j,count)) - 2 / dy^2;
            A(m, getIndex(P,i,j-1,count)) = C(5);
            temp = temp + C(4);
        end
        if (i * dx)^2 + ((j-1) * dy)^2 > 1
            tau = abs(-sqrt(1 - (i*dx)^2) - j*dy) / dy;
            A(m, getIndex(P,i,j,count)) =  A(m, getIndex(P,i,j,count)) - (3 - tau) / (tau * dy^2);
            A(m, getIndex(P,i,j+1,count)) = 2 * (2 - tau) / ((tau + 1) * dy^2);
            A(m, getIndex(P,i,j+2,count)) = (tau - 1) / ((tau + 2) * dy^2);
            temp = temp + 6 / ((tau * (tau + 1) * (tau + 2)) * dy^2);
        elseif (i * dx)^2 + ((j-1) * dy)^2 == 1.0
            A(m, getIndex(P,i,j,count)) = A(m, getIndex(P,i,j,count)) - 2 / dy^2;
            A(m, getIndex(P,i,j+1,count)) = C(4);
            temp = temp + C(5);
        end
    end
    b(m) = temp;
end

u = A \ (f - b);
result = max(abs(u - u_tilde));
s = 'The maximum error in the numerical solution for n = %d is: %1.8f.';
fprintf(s, n, result);