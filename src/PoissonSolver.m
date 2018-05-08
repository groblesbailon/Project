%German Robles
%1456165
%2D Poisson Equation
%May 9th, 2018
%Scientific Computing

clearvars
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

h = (1/(Nx+1))^2;

hx = linspace(ax,bx,Nx);        %equally spaced vector in x direction
hy = linspace(ay,by,Ny);        %equally spaced vector in y direction

[x,y] = meshgrid(hx,hy);        %filled matrix from -pi to pi in both direction to simulate rectangle
y = flipud(y);                  %flipped y vector to have it go from -pi to pi in vertical direction

fa = ((x-ax).^2).*cos(pi.*x/ax);
ga = x.*((x-ax).^2);
%F = cos((pi/2).*(2.*((x-ax)/(bx-ax))+1)).*sin(pi.*((y-ay)/(by-ay)));                    %right hand side
F = zeros(Nx,Ny);
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
    u(i,1) = (1/4)*(2*u(i,1)+u(i-1,1)+u(i+1,1)+(h)*F(i,1));
end

% corners

u(1,1) = (u(1,2)+u(2,1))/2;         %top left
u(1,Ny) = (u(1,Ny-1)+u(2,Ny))/2;    %top right
u(Nx,1) = (u(Nx-1,1)+u(Nx,2))/2;    %bottom left
u(Nx,Ny) = (u(Nx,Ny-1)+u(Nx-1,Ny))/2;

%% gauss seidel 
    
error = 1;
tole = 1e-6;
gaussit = 0;
ukp1 = u;

while error > tole
    gaussit = gaussit + 1;
    for j = 2:Ny-1
        for i = 2:Nx-1
            ukp1(i,j)=.25*(ukp1(i-1,j)+u(i+1,j)+ukp1(i,j-1)+u(i,j+1)+h*F(i,j));
        end
    end
    error =(1/(Nx*Ny))*sum(sum(abs(ukp1-u)));
    u = ukp1;
end

disp('Gauss Iteration for F =')
disp(gaussit) %shows the total number of iteration to converge 

figure(1)
mesh(x,y,u)
xlabel('x')
ylabel('y')
zlabel('U')
title('3D Solution using Gauss Seidel Method for F ')
    
%% sor method

error = 1;
tole = 1e-6;
gaussit = 0;
ukp1 = u;
w = 1.5;

while error > tole
    gaussit = gaussit + 1;
    for j = 2:Ny-1
        for i = 2:Nx-1
            ukp1(i,j)=((w/4)*(u(i+1,j)+ukp1(i-1,j)+ u(i,j+1)+ ukp1(i,j-1)+(h*F(i,j))))+(1-w)*u(i,j);
        end
    end
    error =(1/(Nx*Ny))*sum(sum(abs(ukp1-u)));
    u = ukp1;
end

disp('Successive over Relaxation for F =')
disp(gaussit) %shows the total number of iteration to converge 

figure(2)
mesh(x,y,u)
xlabel('x')
ylabel('y')
zlabel('U')
title('3D Solution using SOR Method for F ')


    