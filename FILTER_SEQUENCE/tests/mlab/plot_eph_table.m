clear all;
close all;


f=load('../check_interpolation.fla');
t1=f(:,1);
equ1 = f(:,2:7);

f=load('../check_nodes.fla');
t2=f(:,1);
equ2 = f(:,2:7);

figure(1)
for i=1:6
subplot(2,3,i)
plot(t1,equ1(:,i),'or')
hold on
plot(t2,equ2(:,i),'.b')
end
