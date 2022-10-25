function [A,x] = LSfit(time, data, f)
%This function will least squares fit a chunk of data using An number of
%sines and cosines
%Inputs: time --> time vector in dnum units of [days]
%        data --> data vector of same length as [time]
%        f --> vector containing frequencies in [1/days] of sine/cosines to fit. The
%        model wil fit length(f) sine/cosine paris with frequencies f 
%Returns: A --> matrix containing the functions to fit to the data
%         x --> vector containing fit coefficients solved for in the least
%         squares method
%
%           To plot data with LSfit: plot(time,data,time,A2*x2)
%written by Charlotte Bellerjeau on 10/12/22

%checks 

[m,n]=size(data);
if m==1 && n>1
    data = data';
elseif m~=1 && n~=1
    print('Argument data must be nx1 column vector')
end

[mm,nn]=size(time);
if mm~=m || nn~=n
    print('Arguments data and time must have matching dimensions')

elseif mm==1 && nn>1
    time = time';
elseif mm~=1 && nn~=1
    print('Argument time must be nx1 column vector')
end

%build A matrix of functions to fit to data
A = ones(length(time),1);
for i=1:length(f)
    A = [A sin(2*pi*time*f(i)) cos(2*pi*time*f(i))];
end

%solve for length(f)+1 coefficients for functions in A
x=(A'*A)\A'*data;

end