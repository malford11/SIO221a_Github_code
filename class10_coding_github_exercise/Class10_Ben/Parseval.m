function [p,variance,sum_spec] = Parseval(DataOG,spectra,df)
%PARSEVAL checks for Parseval's Th. Should equal approximately 1.
variance=nanmean(detrend(DataOG,1).^2);
sum_spec=sum(spectra)*df;
p = sum_spec / variance;
end

