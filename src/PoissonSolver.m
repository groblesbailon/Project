%German Robles
%1456165
%2D Poisson Equation
%May 9th, 2018
%Scientific Computing

clear all
clc
close all

%% variables 

ax = -pi;
ay = -pi;
bx = pi;
by = pi;
Lx = bx - ax;
Ly = by - ay;

Nx = 200;                %nodes in the x direction
Ny = 200;                %nodes in the y direction

hx = linspace(ax,bx,Nx);        %equally spaced vector in x direction
hy = linspace(ay,by,Ny);        %equally spaced vector in y direction

[x,y] = meshgrid(hx,hy);        %filled matrix from -pi to pi in both direction to simulate rectangle
y = flipud(y);                  %flipped y vector to have it go from -pi to pi in vertical direction

fa = ((x-ax).^2).*cos(pi.*x/ax);
ga = x.*((x-ax).^2);
F = cos((pi/2).*(2.*((x-ax)/(bx-ax))+1)).*sin(pi.*((y-ay)/(by-ay)));                    %right hand side
uby = fa;                                                                               %top BC
uay = ga;                                                                               %bottom BC
ubx = (bx.*((bx-ax).^2))+(((y-ay)/(by-ay)).*((((bx-ax).^2).*cos(pi.*bx/ax))-(bx.*((bx-ax).^2))));   %Right BC


