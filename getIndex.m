function r = getIndex(P, i, j, n)
r = 0;
for k = 1:n
    if P(k,1) == i && P(k,2) == j
        r = k;
        break;
    end
end
