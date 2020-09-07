function q=getLandmarksMinMax(data,n,N)
    DSamples=pdist2(data,data,'euclidean');
    q=zeros(1,n);
    q(1) = randperm(N,1);
    qRemaining=setdiff(1:N,q);
    for i=2:n
     [~,indexOfMaxOfMins]=max(min(DSamples(q(1:i-1),qRemaining)));
     q(i)= qRemaining(indexOfMaxOfMins);
    end
end