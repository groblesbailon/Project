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

Nx = 100;                %nodes in the x direction
Ny = 100;                %nodes in the y direction
dx = (Lx)/(Nx+1);
dy = (Ly)/(Ny+1);
h = dx*dy;

hx = linspace(ax,bx,Nx);
hy = linspace(ay,by,Ny);

error = 1;
tolerance = 1e-6;
gaussiter = 0;

u = zeros(Nx,Ny);
ukp1 = u;
F = zeros(Nx-1,Ny-1);
fa = zeros(Nx-1,Ny-1);
ga = zeros(Nx-1,Ny-1);
gb = zeros(Nx-1,Ny-1);

for j=2:Ny-1
    for i=2:Nx-1
        F(i,j) = cos((pi/2)*((2*((hx(i)-ax)/(bx-ax)))+1))*sin(pi*((hy(j)-ay)/(by-ay)));
        fa(i,j) = ((hx(i)-ax)^2)*cos(pi*hx(i)/ax);
        ga(i,j) = hx(i)*((hx(i)-ax)^2);
        gb(i,j) = ga(i,j)+(((hy(j)-ay)/(by-ay))*(fa(i,j)-ga(i,j)));
    end
end

while error > tolerance
    gaussiter = gaussiter + 1;
    for j = 2:1:Ny-1
        for i = 2:1:Nx-1
            ukp1(i,j) = 0.25*(ukp1(i-1,j)+u(i+1,j)+ukp1(i,j-1)+u(i,j+1)-h*F(i,j));
        end
    end
    error=(1/(Nx*Ny))*sum(sum(abs(ukp1-u)));
    u = ukp1;
end
