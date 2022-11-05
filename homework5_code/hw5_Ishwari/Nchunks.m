function [N] = Nchunks(overlap, L, time)
%input L is in n_days in each chunk
%a is overlap fraction
%time is a datetime vector
%output is a vector [nchunks, npoints]
N(1) = (((datenum(time(end)) - datenum(time(1)))/L) - overlap)/(1-overlap);
N(2) = find(datenum(time) == datenum(time(1)) + L);