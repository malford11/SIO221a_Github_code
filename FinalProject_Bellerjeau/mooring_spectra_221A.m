function [mooring] = mooring_spectra_221A(mooring,chunk)


%size of dataset
[m,n] = size(mooring.u);  %[time, depth]

%create up canyon frequency spectrum for each depth using spectrumCB
dtrend = 1;
Hanning = 1;
cosine = 0;
%frk = 1/100;
for i=1:n
    %make spectrum for each depth
    [k,F(:,i),Pars(i)] = spectrumCB(mooring.time,mooring.upSc(:,i),chunk,dtrend,Hanning,cosine);
    %smooth spectrum at each depth
    %Fs(:,i) = specsmoothCB(k,F(:,i),frk);
    Fs = F;
end
mooring.Fs = Fs;
mooring.k = k;
mooring.Pars = Pars;

end