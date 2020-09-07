function dataNormalized=normalizeData(data)
    %Normalize Sample Data
    dataNormalized=zeros(size(data));
    for j=1:size(data,2)
    dataNormalized(:,j)=(data(:,j)-min(data(:,j)))/( max(data(:,j))-min(data(:,j)) );
    end
end