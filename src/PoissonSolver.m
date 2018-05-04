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