function [fa] = alias(fs,fq)
%%% by Charlotte Bellerjeau on 11/25/22 for SIOC 221A
%%% This function computes the aliased frequency of a query frequency
%%% sampled at a given sampling frequency
%%%     Inputs:     fs --> sampling frequency
%%%                 fq --> query frequency
%%%     Output:     fa --> alias frequency


fNyq = fs/2;
M = floor(fq/fNyq);
if ~mod(M,2) %if M is even
    fa = fq - M*fNyq;
else   %if M is odd
    fa = fNyq - (fq - M*fNyq);
end



end