function Alias2(f_signal,dt)
%Alias algorithm. f_signal in hours and dt in days.
f_signal = 1/(f_signal/24); %cycles per day
f_n = 1/2/dt;
M = floor(f_signal/f_n);
if rem(M,2)>0 %if M is odd
    delta = f_signal - M*f_n;
    f_alias = f_n - delta;
elseif rem(M,2)==0 %if M is even
    delta = f_signal - M*f_n;
    f_alias = delta; %cycles per hour
end
f_alias = 1/f_alias; %hours in 1 cycle
% f_alias = f_alias / 24; %days

x = ['The ', num2str(f_signal), ' CPD oscillation sampled at '...
    , num2str(dt), '-day intervals aliases to ', num2str(f_alias)...
    , ' days'];
disp(x)
end

