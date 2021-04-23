function Score = IGD(pbest,NDSet)
% <metric> 

% Copyright 2015-2016 BIMK group

    Distance = zeros(size(NDSet,1),size(pbest,1));
    for i = 1 : size(NDSet,1)
    	Distance(i,:) = sqrt(sum((repmat(NDSet(i,:),size(pbest,1),1)-pbest).^2,2))';
    end
    Distance = min(Distance,[],2);
    Score    = mean(Distance);
end