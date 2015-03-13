clear all;
close all;


addpath(genpath('chebfun'));


f=load('nodes.fla');
t1=f(:,1);
equ1 = f(:,2:7);

f=load('eph_table.fla');
t2=f(:,1);
equ2 = f(:,2:7);

figure(1)
for i=1:6
    subplot(2,3,i)
plot(t1,equ1(:,i),'or')
hold on
plot(t2,equ2(:,i),'.b')
end


x = t1;
a = x(1); b = x(end);

figure(2)

for i=1:6

y = equ1(:,i);
f = polyfit(x,y,100,domain(a,b));

subplot(3,2,i)
plot(x,y,'xk','markersize',12)
hold on
plot(f,'r','linewidth',2)
%title('Discrete polynomial least-squares fit','fontsize',16)

end