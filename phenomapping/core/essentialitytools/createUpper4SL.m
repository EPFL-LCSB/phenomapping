function grRateKOmatrixMirroed = createUpper4SL(grRateKOmatrix, geneList)

grRateKOmatrixMirroed = grRateKOmatrix;

msol = max(max(grRateKOmatrix));
sol = max(max(grRateKOmatrix(grRateKOmatrix<max(max(grRateKOmatrix))-1)));

grRateKOmatrixMirroed = grRateKOmatrix - msol*eye(length(geneList)) + sol*eye(length(geneList));

if grRateKOmatrixMirroed(1,2) == msol
    for i = 1:length(geneList)
        grRateKOmatrixMirroed(1:i,i) = grRateKOmatrixMirroed(i,1:i)';
    end
elseif grRateKOmatrixMirroed(2,1) == msol
    for i = 1:length(geneList)
        grRateKOmatrixMirroed(i,1:i) = grRateKOmatrixMirroed(1:i,i)';
    end
else
    error('no need to mirror the matrix')
end