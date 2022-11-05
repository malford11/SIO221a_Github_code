function [A] = chunkify(overlap, N, data)
%N = [nchunks, npoints] from Nchunks function
%overlap is overlap fraction
%data is data to be chunkified
%output is an npoints*nchunks matrix 

n = floor(N(1));
del = N(2)*(1-overlap);
L = length(data);
A = [];
for j = 1:n
   if  floor((j-1)*del)+ N(2) > L
       A(1: N(2), j) = horzcat(data(1 + floor((j-1)*del): L), zeros(1,floor((j-1)*del)+ N(2) - L))';
   else
       A(1: N(2), j) = data(1 + floor((j-1)*del): floor((j-1)*del)+ N(2))';
   end    
end   