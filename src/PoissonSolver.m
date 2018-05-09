%German Robles
%1456165
%2D Poisson Equation
%APc1-2 version
%May 9th, 2018
%Scientific Computing

clearvars
clc
close all

%% variables that will be used in this problem to solve for 2D poisson and laplace equation

%rectangle of interest
ax = -pi;               %left boundary of rectangle of interest
ay = -pi;               %bottom boundary of rectangle of interest
bx = pi;                %right boundary of rectangle of interest
by = pi;                %top boundary of rectangle of interest
Lx = bx - ax;           %Lx=2pi since thats magnitude of the rectangle in x-direction
Ly = by - ay;           %Ly=2pi since thats magnitude of the rectangle in y-direction

%number of nodes analyzed in this problem
disp('Enter the number of nodes to be analyzed:')       %user friendly input depending on the number of nodes
N = input('Nodes: ');               %choose a desired number of poles

Nx = N;                 %Ny is set equal to the desired number of nodes
Ny = N;                 %Ny is set equal to the desired number of nodes

h2 = (1/(Nx+1))^2;      %since h is the same for deltax and deltay, h^2 is used for the following calculation

hx = linspace(ax,bx,Nx);        %equally spaced vector from -pi to pi using N nodes in x direction
hy = linspace(ay,by,Ny);        %equally spaced vector from -pi to pi usinf N nodes in y direction

[x,y] = meshgrid(hx,hy);        %created matrix from -pi to pi in both direction to simulate rectangle
y = flipud(y);                  %flipped y vector to have it go from -pi to pi in vertical direction

%boundary conditions defined
fa = ((x-ax).^2).*cos(pi.*x/ax);            %equation used for the top boundary condition (dirichlet)
ga = x.*((x-ax).^2);                        %equation used for the bottom boundary condition (dirichlet)
F = cos((pi/2).*(2.*((x-ax)/(bx-ax))+1)).*sin(pi.*((y-ay)/(by-ay)));                    %right hand side for poissons equation (first case) 
%F = zeros(Nx,Ny);              %right hand side for laplaces equation (second case)
uby = fa;                                     %top BC
uay = ga;                                     %bottom BC
ubx = (bx.*((bx-ax).^2))+(((y-ay)/(by-ay)).*((((bx-ax).^2).*cos(pi.*bx/ax))-(bx.*((bx-ax).^2))));   %Right BC

%% boundary conditions prescribed on u matrix

%solvng for outer edges of the u matrix solution
u = zeros(Nx,Ny);       %preallocating the matrix for u matrix solution
u(1,2:Ny-1)=uby(1,2:Ny-1);      %plugging in top BC's defined earlier on the u matrix
u(Nx,2:Ny-1)=uay(Nx,2:Ny-1);    %plugging in bottom BC's defined earlier on the u matrix
u(2:Nx-1,Ny)=ubx(2:Nx-1,Ny);    %plugging in right BC's defined earlier on the u matrix

% neumann boundary condition on the left side of the u matrix

for i = 2:Nx-1
    u(i,1) = (1/4)*(2*u(i,1)+u(i-1,1)+u(i+1,1)+(h2)*F(i,1));    %algorithm used to calculate left BC since a ghost node is requiered 
end

% to provide bettter results, the corners will be calculated as averages of
% the neighboring grids

u(1,1) = (u(1,2)+u(2,1))/2;         %average for top left corner
u(1,Ny) = (u(1,Ny-1)+u(2,Ny))/2;    %average for top right corner
u(Nx,1) = (u(Nx-1,1)+u(Nx,2))/2;    %average for bottom left corner
u(Nx,Ny) = (u(Nx,Ny-1)+u(Nx-1,Ny))/2;   %average for bottom right corner

%% first iterative method used to solve poisson equation: gauss seidel 
    
error = 1;                      %define the error
tole = 1e-6;                    %tolerance used when iterating using the gauss seidel method
gaussit = 0;                    %keeps track of the number of iterations to date
ukp1 = u;                       %defines the u (k+1) matrix to the original u matrix

while error > tole
    gaussit = gaussit + 1;      %iterator count keeps growing for each iteration
    for j = 2:Ny-1
        for i = 2:Nx-1
            ukp1(i,j)=.25*(ukp1(i-1,j)+u(i+1,j)+ukp1(i,j-1)+u(i,j+1)+h2*F(i,j));    %algorithm for the GS method
        end
    end
    error =(1/(Nx*Ny))*sum(sum(abs(ukp1-u)));           %calculated error used for the next iteration
    u = ukp1;
end

%iterations needed to solve
disp('Number of Iterations using the Gauss Seidel method for F=0 =')
disp(gaussit)                   %shows the total number of iteration to converge 

%plots solution
figure(1)                       %plots the u solution on a x, y, and z plane
mesh(x,y,u)
xlabel('x')
ylabel('y')
zlabel('u')
title('Solution using the Gauss Seidel Method for F=0 ')
    
%% second iterative method used to solve poisson equation: SOR method

error = 1;
tole = 1e-6;
gaussit = 0;
ukp1 = u;
w = 1.1;                        %relaxation factor used to solve the SOR method. ideally this value must be within 1 to 2. 
%for this first calculation I used 1.5 but this value needs to be optimzed
%to produce the best results possible in the form of the fastest
%convergence

while error > tole
    gaussit = gaussit + 1;
    for j = 2:Ny-1
        for i = 2:Nx-1
            ukp1(i,j)=((w/4)*(u(i+1,j)+ukp1(i-1,j)+ u(i,j+1)+ ukp1(i,j-1)+(h2*F(i,j))))+(1-w)*u(i,j);   %algorithm for the SOR method
        end
    end
    error =(1/(Nx*Ny))*sum(sum(abs(ukp1-u)));
    u = ukp1;
end

%iterations needed to solve
disp('Number of Iterations using the SOR method for F=0 =')
disp(gaussit) %shows the total number of iteration to converge 

%plots solution
figure(2)                           %plots a second figure to display the SOR solution on a x, y, and z plane
mesh(x,y,u)
xlabel('x')
ylabel('y')
zlabel('u')
title('Solution using the SOR Method for F=0 ')


    