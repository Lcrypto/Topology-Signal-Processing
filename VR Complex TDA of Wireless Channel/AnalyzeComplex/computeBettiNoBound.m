function betti=computeBettiNoBound(delK,delKPlusOne)
m=1;

delKReduced=SNF(delK);

kernelDimension=size(delK,2)-rank(delK);
betti=kernelDimension;

end