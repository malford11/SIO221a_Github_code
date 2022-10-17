function statistics_variance = stats_variance(var1,var2)
% This function will calculate the variance of two datasets with standard
% errors

N1 = length(var1);
N2 = length(var2);
Var1 = std(var1)^2;
Var2 = std(var2)^2;
Varerr1=Var1*sqrt(2/(N1-1));
Varerr2=Var2*sqrt(2/(N2-1));

statistics_variance=[Var1, Var2, Varerr1, Varerr2];
end