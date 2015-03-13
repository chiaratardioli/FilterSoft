clear all; close all;

addpath(genpath('chebfun'));


x = chebfun('x');
format long
disp('Cheb coeffs of 99x^2 + x^3:')
p = 99*x.^2 + x.^3;
a = chebcoeffs(p)

disp('Cheb coeffs of exp(x):')
a = chebcoeffs(exp(x))

FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth';
plotcoeffs(exp(x),'.-',LW,1,MS,20), grid on
xlabel('degree n',FS,14)
ylabel('|a_n|',FS,14), ylim([1e-17 1e1])
title('Chebyshev coefficients of exp(x)',FS,14)

plotcoeffs(exp(x)./(1+10000*x.^2)), grid on
xlabel('degree n',FS,12), ylabel('|a_n|',FS,12)
ylim([1e-18 1])
title('Chebyshev coefficients of exp(x)/(1+10000x^2)',FS,14)

f = sign(x);
figure, plot(f,'k',LW,2), ylim([-1.5 1.5])
title('sign(x)',FS,14)

p = chebfun(f,'trunc',10);
a = chebcoeffs(p)

hold on
plot(p,'m',LW,2)
title('sign(x) and truncated Chebyshev series',FS,14)

pinterp = chebfun(f,10);
plot(pinterp,'--','color',[0 .8 0],LW,2)
title('Same, also with Chebyshev interpolant',FS,14)

npts = 100;
x = linspace(-1,1,npts);
y = 1./(1+25*x.^2) + 1e-1*randn(1,npts);
f = polyfit(x,y,10,domain(-1,1));
plot(x,y,'xk','markersize',12)
hold on, plot(f,'r','linewidth',2)
title('Discrete polynomial least-squares fit','fontsize',16)