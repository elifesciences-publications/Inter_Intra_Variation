function [norm_X] = autoscale(X)
%[norm_X] = autoscale(X)
%   Detailed explanation goes here

for D = 1:size(X,2)
    SD_X = std(X(:,D));
    mu_X = mean(X(:,D));
    for obs = 1:size(X,1)
        norm_X(obs,D) = (X(obs,D)-mu_X)/SD_X;
    end
end
end

