function [x, errorbar] = ErrorBar(frq,f_vector,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Confidence Bounds
nu=10;
err_low = nu/chi2inv(1-.05/2,nu);
err_high = nu/chi2inv(.05/2,nu);

m = find(f_vector == frq);
x = ones(10,1) * f_vector(m);
errorbar = linspace(err_low*y(m),err_high*y(m),10);
end

