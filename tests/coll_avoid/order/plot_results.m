clear all
close all
clc

re = 6.378137e+003; % [km]
UD=42164.1697748545;

f=load('results_cart.dat');
t=f(:,1)*2*pi;
x1=f(:,2)*UD;
x2=f(:,3)*UD;
x3=f(:,4)*UD;


figure(1)

[xe,ye,ze] = sphere;
plot3(xe*re,ye*re,ze*re,'b')
hold on
plot3(x1,x2,x3,'.')