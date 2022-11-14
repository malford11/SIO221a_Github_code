function [Data] = CreateData(n,m,z,dt)
%Create Data with given number of elements for Part 1
%Normalize data:
% "so we could divide our spectrum by twice the Nyquist frequency to have
% energy in units appropriate for comparing if we wanted to have our
% integrals match." - Lecture 10
Data = struct;
Data.WhiteNoise = randn(n,m,z); %Gaussian white noise dataset
Data.TimeInt = cumsum(Data.WhiteNoise,1) * dt; %Time Integral
                                                %dt=10 min
                                               

Data.Names = {'White Noise', 'Time Integral'};
end

