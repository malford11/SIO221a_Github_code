function statistics_mean = stats_mean(var1,var2)
% This function will calculate the mean of two datasets with standard
% errors

N1 = length(var1);
N2 = length(var2);
Mean1= mean(var1);
Mean2= mean(var2);
Meanerr1=std(var1)/sqrt(N1);
Meanerr2=std(var2)/sqrt(N2);

statistics_mean=[Mean1, Mean2, Meanerr1, Meanerr2];
end