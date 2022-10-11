function h=freqline(fi,sty)
%function freqline(fi,sty)
%Place a vertical line at the location specified by fi.
%
%MHA 03/04
%01/05: return the handle.
%
if nargin < 2
sty='k--';
end

xv=[fi fi]';
yl=ones(size(xv))*ylim;
hold on
h=loglog(xv,yl',sty);
hold off
