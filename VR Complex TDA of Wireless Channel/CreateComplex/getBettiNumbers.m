function betti=getBettiWitness(data,landmarks,v,simplexDimension,R,plotOn)
dim=size(data,2);

D=pdist2(data(landmarks,:),data,'euclidean');

%1) Build Simplicial Complex 
[graph,complex]=getGraphAndComplex(D,R,v,simplexDimension);

%Plotting
if plotOn==1 
    figure
    for k=1:size(complex,2)
        fprintf('#%1d-Simplicies: %1d\n',k-1,size(complex{k},1));
    end
    hold on
    gplot(graph,data(landmarks,:),'ro-') %Plot NeighbourhoodGraph
       if dim==3
        plot3(data(:,1),data(:,2),data(:,3),'o','Color','b');
    elseif dim<=2
        plot(data(:,1),data(:,2),'o','Color','b')
    end
end

%2) Create a "Boundary Matrix" for each order k.
boundaryMatricies=computeBoundaryMatricies(complex,k);
set(0,'RecursionLimit',1500)

%3) Get the Betti Numbers by Analyzing the Boundary Matrices.
if size(complex,2)>1
for k=1:size(complex,2)-1
betti(k)=computeBetti(boundaryMatricies{k},boundaryMatricies{k+1});
end
betti(k+1)=computeBettiNoBound(boundaryMatricies{k+1});
else
betti(1)=computeBettiNoBound(boundaryMatricies{1});
end
end