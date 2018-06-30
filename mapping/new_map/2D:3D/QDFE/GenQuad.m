function quad = GenQuad(order)

% Load and read QDFE quadrature
fid     = fopen('QDFEQUAD.DAT');
quadAll = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f');
fclose(fid);

% First line of quadrature data
lineStart = 0; 
for i = 0 : order - 1
    orderUse = 3 ^ i - 1;
    numSubSq = (orderUse + 1) ^ 2 * 9 * 3;
    lineStart = lineStart + numSubSq;
end
lineStart = lineStart + 1;

% Number of quadrature directions 
orderUse = 3 ^ order - 1;
numDirs = (orderUse + 1) ^ 2 * 9 * 3;

% Resize storage
quad = cell(numDirs, 1);

% Store quadrature
for i = 1 : numDirs
    for j = 1 : 13
        quad{i}(j) = quadAll{j}(lineStart + i - 1);
    end
end

end