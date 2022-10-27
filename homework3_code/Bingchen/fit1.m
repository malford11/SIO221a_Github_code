function [data_fit,x] = fit1(data,num_date,P)

A_fit=[ones(length(num_date),1) sin(2*pi*num_date/P) ...
    cos(2*pi*num_date/P)];
x = inv(A_fit'*A_fit)*A_fit'*data;
data_fit = A_fit*x;
x(4) = norm(x(2:3));

end 