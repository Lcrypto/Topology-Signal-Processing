function lowerNeighbourIndicies=lowerNeighbours(G,u)
    lowerNeighbourIndicies=find(G(1:u,u)==1);
end