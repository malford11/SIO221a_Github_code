function [stats] = stats(data,var)
N = length(data.(var));
mean = nanmean(data.(var));
std = nanstd(data.(var));
std_error = std/sqrt(N);
variance = std^2;
std_var = (variance * sqrt(2/(N-1)));

stats = [mean, std, std_error, variance, std_var];
end

