function [TL] = SpectraPlots(n,m)
%Create nxm tiled layout graph axes with labels
TL = tiledlayout(n,m);
xlabel(TL,'Cycles per Day');
ylabel(TL,'(measurement units)^2/ (cpd) normalized)');
end

