function [f, spectrum, p] = spectrumIsh(time, data, n, overlap, flag)
%function for computing the spectrum from a given data series
%with options for overlap and Hann window
%ishwari %23 Nov 22

%%% inputs %%%
%time in number of days since..
%n = number of segments 
%overlap = overlap fraction %ex: for 50% overlap, overlap = 0.5 
%flag == 1, hanning window
%flag == 0, no window
%no input for detrending, it WILL happen

%%% outputs %%%
%spectrum = spectrum
%f = frequency vector
%p = ratio test for Parseval

%% chunkifying and windowing data

tdiff = nanmean(diff(time));        %time difference in days between 2 data points
L = length(data);                   %total number of points
l = L/(n*(1-overlap) + overlap);    %number of points in each chunk
l = floor(l);

%loop for chunkifying data 
del = l*(1-overlap);
A = [];
for j = 1:n
   if  floor((j-1)*del)+ l > L
       A(1: l, j) = horzcat(data(1 + floor((j-1)*del): L), zeros(1,floor((j-1)*del)+ l - L))';
   else
       A(1: l, j) = data(1 + floor((j-1)*del): floor((j-1)*del)+ l)';
   end    
end  
%matrix A is an l*n matrix = length of chunk * number of chunks matrix


%loop for chunkifying time
T = [];
for j = 1:n
   if  floor((j-1)*del)+ l > L
       T(1: l, j) = horzcat(time(1 + floor((j-1)*del): L), zeros(1,floor((j-1)*del)+ l - L))';
   else
       T(1: l, j) = time(1 + floor((j-1)*del): floor((j-1)*del)+ l)';
   end    
end  

%windowing
if flag == 1
    hanwin = cos(pi*T/ (l*tdiff)).^2;             %defining Hanning window
    hanwin = hanwin/sqrt(8/3);                    %normalizing
    A(1:end,:) = hanwin(1:end,:).*A(1:end,:);     %convolving
end
A = detrend(A);

%% defining frequency vector

df = 1/(l*tdiff);           %fundamental frequency, depends on chunk length
fn = 1/(2*tdiff);           %Nyquist freq    
f = 0:df:fn;                %freq vector
%all in cpd

%% computing spectrum

Z = fft(A);
spectrum = abs(Z(1:floor(l/2)+1,:)).^2; %throwing away bottom half, squaring
spectrum = (2/l^2))*spectrum;           %factor of 2, dividing by N^2
spectrum = spectrum/df;                 %defining spectrum

spectrum = nanmean(spectrum,2);         %averaging over chunks

%% parseval

specSum = sum(spectrum);
var = std(detrend(data))^2;
p = specSum/var;


end
