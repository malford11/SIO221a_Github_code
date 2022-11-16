t=0:.2:30;
x=0:.2:300;

[tt,xx]=meshgrid(t,x);

k=2*pi/150;
g=9.8;
om=-sqrt(g*k);

eta=sin(k*xx - om*tt);

%mesh(t,x,eta)

for c=2:length(t)
    
    subplot(211)
    plot(x,eta(:,c))
    xlim([0 max(x)])
    subplot(212)
    ezpc(t(1:c),x,eta(:,1:c));axis xy
    xlim([0 max(t)])
    ylabel('x / m')
    xlabel('time / s')
    drawnow
    shg
end
