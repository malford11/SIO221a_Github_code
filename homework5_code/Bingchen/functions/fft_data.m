function [f,xs,spec]=fft_data(data,dt)
    
    N = length(data);
    T=dt*N;
    df = 1/T;
    fn = 1/2/dt;
    f=-fn:df:fn-df; 

    x=fft(data);
    x_norm =fft(data)/sqrt(N);
    xs= fftshift(x_norm);

    spec = (xs.*conj(xs)).*2/df; %Since we only plot half of the spectrum, we correct for the lost variance by *2.

end 