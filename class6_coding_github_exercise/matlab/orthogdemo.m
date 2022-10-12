function [out] = orthogdemo(n,m)
% inputs: n --> frequency of 1st sine wave
%        m --> frequency of 2nd sine wave
% returns: out struct which contains
%           time --> time vector with same units as freq
%           wave1 --> sine wave with frquency n
%           wave2 --> sine wave with frequency m
%           product --> wave1*wave2
%           integral --> integral of product over time vector

    out.time = 0:0.01:1;
    out.wave1 = sin(2*pi*n.*out.time);
    out.wave2 = sin(2*pi*m.*out.time);
    out.product = out.wave1.*out.wave2;
    out.integral = trapz(out.time,out.product);
    out.n;
    out.m;
end