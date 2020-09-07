function betti=computeBetti(delK,delKPlusOne)
betti=1;

delKReduced=SNF(delK);

kernelDimension=size(delKReduced,2)-rank(delKReduced);


delKPlusOneReduced=SNF(delKPlusOne);
imageDimension=rank(delKPlusOneReduced);

betti=kernelDimension-imageDimension;



end