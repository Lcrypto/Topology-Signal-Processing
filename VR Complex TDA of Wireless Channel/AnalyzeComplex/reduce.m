function del=reduce(del,x)
% Algorithm from "Computational Topology: An Introduction, 2008 (H.Edelsbruneer &
% J.Harer)"
[r,c]=find(del==1);
rcCombosGreaterThan_x=find(((r>=x) & (c>=x)));
if ~isempty(rcCombosGreaterThan_x)
    k=r(rcCombosGreaterThan_x(1));
    l=c(rcCombosGreaterThan_x(1));

   %Swap Rows
   del([x k],:) = del([k x],:);

   %Swap Columns
   del(:,[l x]) = del(:,[x l]);
   
    for i=x+1:size(del,1)
        if (del(i,x)==1)
           del(i,:)=del(x,:)+del(i,:) ;
        end
    end

    for j=x+1:size(del,2)
        if (del(x,j)==1)
           del(:,j)=del(:,x)+del(:,j) ;
        end
    end
    del=reduce(del,x+1);
end

return
end