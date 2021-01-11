m = 10;
% Boundary pts.
b = zeros(m+2, 1);
% Length of evenly spaced mesh.
h = 1 / m;

c = (1/2 * (exp(1) - exp(-1)) - pi);

% % Gaussian quadrature
w = zeros(5,1); % weights
w(1) = 128/225;
w(2) = (322 + 13*sqrt(70)) / 900;
w(3) = (322 + 13*sqrt(70)) / 900;
w(4) = (322 - 13*sqrt(70)) / 900;
w(5) = (322 - 13*sqrt(70)) / 900;
z = zeros(5,1); % points
z(1) = 0;
z(2) = 1/3*sqrt(5 - 2*sqrt(10/7));
z(3) = -1/3*sqrt(5 - 2*sqrt(10/7));
z(4) = 1/3*sqrt(5 + 2*sqrt(10/7));
z(5) = -1/3*sqrt(5 + 2*sqrt(10/7));

int11 = zeros(m+2, 1);
int12 = zeros(m+2, 1);
int13 = zeros(m+2, 1);
int14 = zeros(m+2, 1);

int21 = zeros(m+2, 1);
int22 = zeros(m+2, 1);
int23 = zeros(m+2, 1);
int24 = zeros(m+2, 1);

int31 = zeros(m+1, 1);
int32 = zeros(m+1, 1);
int33 = zeros(m+1, 1);
int34 = zeros(m+1, 1);
for i = -2:m-1
    for k = 1:5
        if i == -2
            x = h/2*(z(k) + 1); % input for x_{i} to x_{i+1}.
            int13(1) = int13(1) + 1/2 * w(k) ...
                * ((pi^2 + 1)*sin(pi*x)*(3*x/h - 9/2*(x/h)^2 ...
                + 7/4*(x/h)^3));
            int23(1) = int23(1) + 1/2 * w(k) ...
                * ((1 + c*x)*(3*x/h - 9/2*(x/h)^2 ...
                + 7/4*(x/h)^3));
            int33(1) = int33(1) + 1/2 * w(k) ...
                * (c*(3/h - 9/h*(x/h) ...
                + 21/(4*h)*(x/h)^2));
            x = h/2*(z(k) + (2*i+3)); % input for x_{i+1} to x_{i+2}.
            int14(1) = int14(1) + 1/2 * w(k) ...
                * ((pi^2 + 1)*sin(pi*x)*(1/4*((2*h-x)/h)^3) ...
                - 1/2*((x - (i)*h)/h)^3);
            int24(1) = int24(1) + 1/2 * w(k) ...
                * ((1 + c*x)*(1/4*((2*h-x)/h)^3);
            int34(1) = int34(1) + 1/2 * w(k) ...
                * (c*(-3/(4*h)*((2*h-x)/h)^2));
        elseif i == -1
            % Implement B_-1 cubic for b vector.
        else
            if i < m-1
                x = h/2*(z(k) + (2*i+1)); % input for x_{i} to x_{i+1}.
                int11(i+1) = int11(i+1) + 1/2 * w(k) * ...
                    ((pi^2 + 1)*sin(pi*x)*((x - (i)*h)/h)^3/6);
                int21(i+1) = int21(i+1) + 1/2 * w(k) * ...
                    ((1 + c*x)*((x - (i)*h)/h)^3/6);
                int31(i+1) = int31(i+1) + 1/2 * w(k) * ...
                    (c*((x - (i)*h)/h)^2/(2*h));
            end
            if i+1 < m-2
                x = h/2*(z(k) + (2*i+3)); % input for x_{i+1} to x_{i+2}.
                int12(i+1) = int12(i+1) + 1/2 * w(k) ...
                    * ((pi^2 + 1)*sin(pi*x)*(2/3 - ...
                    2*(x - (i)*h)/h + 2*((x - (i)*h)/h)^2) ...
                    - 1/2*((x - (i)*h)/h)^3);
                int22(i+1) = int22(i+1) + 1/2 * w(k) ...
                    * ((1 + c*x)*(2/3 - ...
                    2*(x - (i)*h)/h + 2*((x - (i)*h)/h)^2) ...
                    - 1/2*((x - (i)*h)/h)^3);
                int32(i+1) = int32(i+1) + 1/2 * w(k) ...
                    * (c*(-2/h + 4/h*((x - (i)*h)/h) ...
                    - 3/(2*h)*((x - (i)*h)/h)^2));
            end
            if i+2 < m-2
                x = h/2*(z(k) + (2*i+5)); % input for x_{i+2} to x_{i+3}.
                int13(i+1) = int13(i+1) + 1/2 * w(k) ...
                    * ((pi^2 + 1)*sin(pi*x)*(-22/3 + ...
                    10*(x - (i)*h)/h - 4*((x - (i)*h)/h)^2) ...
                    + 1/2*((x - (i)*h)/h)^3);
                int23(i+1) = int23(i+1) + 1/2 * w(k) ...
                    * ((1 + c*x)*(-22/3 + ...
                    10*(x - (i)*h)/h - 4*((x - (i)*h)/h)^2) ...
                    + 1/2*((x - (i)*h)/h)^3);
                int33(i+1) = int33(i+1) + 1/2 * w(k) ...
                    * (c*(10/h - 8/h*((x - (i)*h)/h) ...
                    + 3/(2*h)*((x - (i)*h)/h)^2));
            end
            if i+3 < m-3
                x = h/2*(z(k) + (2*i+7)); % input for x_{i+3} to x_{i+4}.
                int14(i+1) = int14(i+1) + 1/2 * w(k) ...
                    * ((pi^2 + 1)*sin(pi*x)*(1/6*(((i+4)*h - x)/h)^3));
                int24(i+1) = int24(i+1) + 1/2 * w(k) ...
                    * ((1+c*x)*(1/6*(((i+4)*h - x)/h)^3));
                int34(i+1) = int34(i+1) + 1/2 * w(k) ...
                    * (c*(-1/(2*h)*(((i+4)*h - x)/h)^2));
            end
        end
    end
end

% y evaluate at each internal grid pt.
y = zeros(m, 1);

for j = 1:m
    x = j * h;
    y(j) = 1 / 2 * (exp(x) + exp(-x)) + sin(pi * x);
end
% The coefficients of a_ij loaded in by row.
A = sparse(m+2, m+2);

for i = 1:m+2
    for j = 1:m+2
        if i - j == 0
            A(i,j) = 2 / h + 2 * h / 3;
        elseif i - j == 1
            A(i,j) = -1 / h + h * (1 / 6);
        elseif i - j == -1
            A(i,j) = -1 / h + h * (1 / 6);
        end
    end
    b(i) = int11(i) + int12(i) + int13(i) + int14(i) - ...
        (int21(i) + int22(i) + int23(i) + int24(i)) + ...
        int31(i) + int32(i) + int33(i) + int34(i);
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