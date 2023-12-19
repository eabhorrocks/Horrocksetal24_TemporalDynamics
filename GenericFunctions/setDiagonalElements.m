function M = setDiagonalElements(M,val)

M(find(eye(size(M)))) = val;