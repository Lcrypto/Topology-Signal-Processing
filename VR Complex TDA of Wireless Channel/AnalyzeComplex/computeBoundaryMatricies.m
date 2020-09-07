function boundaryMatricies=computeBoundaryMatricies(kSkeletonOfVRComplex,simplexDimension)

for k=1:simplexDimension
    boundaryMatricies{k}=computeBoundaryMatrix(kSkeletonOfVRComplex,k);
end

end



