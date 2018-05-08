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

%% boundary conditions on u matrix

u = zeros(Nx,Ny);       %preallocating for u matrix solution
u(1,2:Ny-1)=uby(1,2:Ny-1); %top BC on u matrix
u(Nx,2:Ny-1)=uay(Nx,2:Ny-1);    %bottom BC on u matrix
u(2:Nx-1,Ny)=ubx(2:Nx-1,Ny);    %right BC on u matrix

% neumann boundary condition

%u(2:Nx-1,1:Ny-1)=0.25*(2*u(2:Nx-1,1:Nx-1)+u(
for i = 2:Nx-1
    u(i,1) = (1/4)*(2*u(i,1)+u(i-1,1)+u(i+1,1)+((1/Nx))*F(i,1));
end

% corners

u(1,1) = (u(1,2)+u(2,1))/2;         %top left
u(1,Ny) = (u(1,Ny-1)+u(2,Ny))/2;    %top right
u(Nx,1) = (u(Nx-1,1)+u(Nx,2))/2;    %bottom left
u(Nx,Ny) = (u(Nx,Ny-1)+u(Nx-1,Ny))/2;

%% gauss seidel 
    
    
    
    
    
    