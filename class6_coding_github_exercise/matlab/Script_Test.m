out = orthogdemo(.5,.5);
figure(1)
clf
tiledlayout(1,3);
%title(m, ['n = ', num2str(out.n), 'm = ', num2str(out.m), 'Integral = ', num2str(out.integral)])
nexttile
plot(out.time,out.wave1)
title('wave 1')
nexttile
plot(out.time,out.wave2)
title('wave 2')
nexttile
plot(out.time,out.product)
title('wave 3')
