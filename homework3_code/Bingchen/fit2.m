function [data_fit,x] = fit2(data,num_date,P)

A_fit=[ones(length(num_date),1) sin(2*pi*num_date/P) ...
    cos(2*pi*num_date/P) sin(4*pi*num_date/P) cos(4*pi*num_date/P)];
x = inv(A_fit'*A_fit)*A_fit'*data;
data_fit = A_fit*x;
x(6) = norm(x(2:5));

end 