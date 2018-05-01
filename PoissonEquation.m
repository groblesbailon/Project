%German Robles
%1456165
%Poisson Equation Project
%Scientific Computing for Mechanical Engineers

clearvars
clc
close all
%% Variables

ax = -pi;
ay = -pi;
bx = pi;
by = pi;
Lx = bx - ax;
Ly = by - ay;

Nx = 50;    %nodes in the x direction
Ny = 50;    %nodes in the y direction
dx = (Lx)/(Nx+1);
dy = (Ly)/(Ny+1);


