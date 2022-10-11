function statistics = stats_CB(data)
% This function will calculate the mean and standard deviation of a dataset
% with error bars
% Input: a vector containing a single data type
% Output: [Mean, Meanerr, Var, Varerr], where
%       Mean = sample mean of data
%       Meanerr = standard error of the mean of data
%       Var = variance of data 
%       Varerr = stadard error of the variance of data
% by Charlotte Bellerjeau on 10/7/2022

N = length(data);
Mean = mean(data);
Var = std(data)^2;
Meanerr=std(data)/sqrt(N);
Varerr=Var*sqrt(2/(N-1));

statistics=[Mean, Meanerr, Var, Varerr];
end
