function [f,spec]=fft_data(data,dt)
    data = detrend(data);
    N = length(data);
    T=dt*N;
    df = 1/T;
    fn = (1/2)/dt;
    %f=-fn:df:fn-df; 
    f = 0:df:fn;                
%% method 1 
%     x=fft(data);
%     x_norm =fft(data)/N;
%     %xs= fftshift(x_norm);
% 
%     xs_pos = x_norm(1:floor(N/2)+1);
%     spec = (xs_pos.*conj(xs_pos)).*2/df;
% 
%     %spec = (xs.*conj(xs)).*2/df; %Since we only plot half of the spectrum, we correct for the lost variance by *2.
%     %spec = (abs(xs).^2)*2/df; %Since we only plot half of the spectrum, we correct for the lost variance by *2.
%     %spec = spec(1:(floor(length(spec)/2)+1));

%% method 2

    a=fft(data);
    a_pos = a(1:floor(N/2)+1);
    amp = (abs(a_pos).^2).*2 / N.^2; 
    spec = amp/df;


%% method 3 from Ben 
%     a=fft(data);
%     amp=abs(a(1:N/2+1)).^2; 
% %     dt=(dt/60)/24;
% %     T=dt*N;
% %     df=1/T;
% %     fn=(1/2)/dt;
% %     f=0:df:fn; %frequency vector, cpd, goes from 0 to Nyquist.
%     amp = amp / N.^2; % first correct for the MATLAB normalization
%     amp = amp .* 2; %we threw out half of the spectrum; so correct for the lost variance.
%     spec = amp / df; % this is then the definition of the spectrum
end 