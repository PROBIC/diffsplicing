function [P R AP] = getPrecRec(values,I)

[sorted_values,rank]=sort(values,'descend');
P=cumsum(I(rank))./(cumsum(~I(rank))+cumsum(I(rank)));  % precision vector
R=cumsum(I(rank))/sum(I);   % recall vector
AP = ((P)'*I(rank))/sum(I); % average precision

end


