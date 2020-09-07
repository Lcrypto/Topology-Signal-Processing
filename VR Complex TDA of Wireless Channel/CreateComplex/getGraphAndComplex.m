function [graph,complex,kMax]=getGraphAndComplex(D,R,v,k)
tic
N=size(D,2);
n=size(D,1);
graph=zeros(n,n);
m=zeros(1,N);

if v>0
   for i=1:N
     column=D(:,i);
     [~,index]=sort(column);
     m(i)=column(index(v)); %vth smallest entry of i column.
   end
end

%1) Get Edges
edges=zeros(0,0);
for ve=1:n
   for l=1:n
       if (l~=ve)
          edge=sort([l ve]);
          if ~(ismember(edge,edges,'rows','legacy'))
              witnessExists=0;
              for i=1:N
                if(max(D(edge(1),i),D(edge(2),i))<=(R+m(i)))
                   witnessExists=1; 
                end
              end
              if (witnessExists==1) 
                edges(end+1,:)=edge;
                graph(l,ve)=1;
                graph(ve,l)=1;
              end
          end
       end
   end
end

%2) Build Higher Simplicies from edges (Using Fast Rips Implementation from Zomorodian).
[complex,kMax]=computeVRComplex(graph,k);
end