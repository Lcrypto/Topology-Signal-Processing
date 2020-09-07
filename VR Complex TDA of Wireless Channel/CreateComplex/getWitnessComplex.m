function [graph,complex]=getWitnessComplex(data,landmarks,v,simplexDimension,R,plotOn)
dim=size(data,2);

D=pdist2(data(landmarks,:),data,'euclidean');
[graph,complex]=getGraphAndComplex(D,R,v,simplexDimension);

if plotOn==1 
    figure
    for k=1:size(complex,2)
        fprintf('#%1d-Simplicies: %1d\n',k-1,size(complex{k},1));
    end
    hold on
    % Plotting
    gplot(graph,data(landmarks,:),'ro-') %Plot NeighbourhoodGraph
        if dim==3
        plot3(data(:,1),data(:,2),data(:,3),'o','Color','b');
    elseif dim<=2
        plot(data(:,1),data(:,2),'o','Color','b')
    end
end

end