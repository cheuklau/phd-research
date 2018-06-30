function surfaceArea = calcSA(xMin, xMax, yMin, yMax)

N = 16;
surfaceArea = 0;
[fx, wx] = lgwt(N, xMin, xMax);
[fy, wy] = lgwt(N, yMin, yMax);
for i = 1 : N
    for j = 1 : N
        surfaceArea = surfaceArea + ...
            (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2)))...
            * wx(i) * wy(j);
    end
end

end