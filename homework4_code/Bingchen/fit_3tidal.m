function [data_fit,x] = fit_3tidal(data,time)

P1 = 25.82;
P2 = 23.93;
P3 = 12.42;

A_fit=[ones(length(time),1) sin(2*pi*time/P1) ...
    cos(2*pi*time/P1) sin(2*pi*time/P2) cos(2*pi*time/P2) ...
    sin(2*pi*time/P3) cos(2*pi*time/P3)];
x = inv(A_fit'*A_fit)*A_fit'*data;
data_fit = A_fit*x;
x(8) = norm(x(2:7));

end 