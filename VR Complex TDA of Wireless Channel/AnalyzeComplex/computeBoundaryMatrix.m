function boundaryMatrix=computeBoundaryMatrix(kSkeletonOfVRComplex,simplexDimension)

KSimplicies=kSkeletonOfVRComplex{simplexDimension};
if (simplexDimension==1)
    boundaryMatrix=zeros(size(KSimplicies,1),1)';
    return;
end


KMinusOneSimplicies=kSkeletonOfVRComplex{simplexDimension-1};
nKSimplicies=size(KSimplicies,1);
nKMinusOneSimplicies=size(KMinusOneSimplicies,1);
boundaryMatrix=zeros(nKMinusOneSimplicies,nKSimplicies);

for i=1:nKMinusOneSimplicies
   for j=1:nKSimplicies
       
      if KSimplicies(j,1)>KMinusOneSimplicies(i,1)
          break
      end
      
      if KMinusOneSimplicies(i,simplexDimension-1)>KSimplicies(j,simplexDimension)
        continue 
      end
     
       if ismembc(KMinusOneSimplicies(i,:),KSimplicies(j,:))
          [a, indiciesOfMissingElements] = find(ismembc(KSimplicies(j,:), KMinusOneSimplicies(i,:))==0);
          boundaryMatrix(i,j)=(-1)^mod(indiciesOfMissingElements+1,2);
       end
   end
end


end
