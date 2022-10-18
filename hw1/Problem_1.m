clear all;
clc;
close all;



%% Known variables of the problem 


% Number of vertices
N = 3;
% Time step size
dt = 0.01; % second
% Rod length
RodLength = 0.1; % meter
% Discrete length
deltaL = RodLength / (N-1);
% Radius of spheres
R1 = 0.005;R2 = 0.005;R3 = 0.005;%metres 
% Density
rho_metal = 7000;% kg/m^3
rho_f = 1000;% kg/m^3
rho = rho_metal - rho_f;% kg/m^3
% Rod radius
r0 = 0.001;% metre
% Young's modulus
Y = 1e9; % Pascals,Using Y instead of E to avoid ambiguity
% Gravity
g = 9.8; % m/s^2
% Viscosity
visc = 1000; % Pa-s
% Total time
totalTime = 10.00; % seconds
% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4;% bending stifness
EA = Y * pi * r0^2;% stretching stifness 
node_number=[1 4 9 99 999];
counter=1;
%% Geometry of the object 


nodes = zeros(N, 2);
for c = 1:N
nodes(c,1) = (c-1) * deltaL; % x-co-ordinates change 
nodes(c,2) = 0; % y co-ordinates are always zero 
end



%% Mass matrix
M = zeros(2*N,2*N);
M(1,1) = 4/3*pi*R1^3*rho_metal;% x1
M(2,2) = 4/3*pi*R1^3*rho_metal;% y1
M(3,3) = 4/3*pi*R2^3*rho_metal;% x2 
M(4,4) = 4/3*pi*R2^3*rho_metal;% y2
M(5,5) = 4/3*pi*R3^3*rho_metal;% x3
M(6,6) = 4/3*pi*R3^3*rho_metal;% y3
%% Viscous damping matrix
C = zeros(6,6);
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;% x1
C(2,2) = C1;% y1
C(3,3) = C2;% x2
C(4,4) = C2;% y2
C(5,5) = C3;% x3
C(6,6) = C3;% y3

%% Gravity/Weight Vector 
W = zeros(2*N,1);
W(2) = -4/3*pi*R1^3*rho*g;% y1
W(4) = -4/3*pi*R2^3*rho*g;% y2
W(6) = -4/3*pi*R3^3*rho*g;% y3

%% Initial DOF vector
q0 = zeros(2*N,1); % old DOFs 
for c=1:N
q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
q0 ( 2*c ) = nodes(c,2); % y coordinate
end
% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector 


%% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros( Nsteps, 1); % y-position of R2
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q(4);
all_mid_v(1) = u(4);
% Tolerance
tol = EI / RodLength^2 * 1e-3;

for c=1:Nsteps-1
%fprintf('Old Time = %f,New Time = %f\n', (c-1) * dt,c*dt );
q = q0; % Guess of the Position
if c==1
    figure(counter)
xcoords=q(1:2:end) %x1,x2,x3
ycoords=q(2:2:end) %y1,y2,y3
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal;
xlabel('Position, m [meter]');
ylabel('Position, m [meter]');
end
% Newton Raphson
err = 10 * tol;
while err > tol
% Inertia
f = M / dt * ( (q-q0) / dt - u );
J = M / dt^2;

% Elastic forces

% Linear spring 1 between nodes 1 and 2
xk = q(1);% x1
yk = q(2);% y1
xkp1 = q(3);% x2
ykp1 = q(4);% y2
dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
f(1:4) = f(1:4) + dF;
J(1:4,1:4) = J(1:4,1:4) + dJ;% Jacobian is not a vector it is a matrix
% Linear spring 2 between nodes 2 and 3
xk = q(3);% x2
yk = q(4);% y2
xkp1 = q(5);% x3
ykp1 = q(6);% y3
dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
f(3:6) = f(3:6) + dF;
J(3:6,3:6) = J(3:6,3:6) + dJ;
% Bending spring between nodes 1, 2, and 3(located at node 2)
xkm1 = q(1);% x1
ykm1 = q(2);% y1
xk = q(3);% x2
yk = q(4);% y2
xkp1 = q(5);% x3
ykp1 = q(6);% y3
curvature0 = 0;
dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
curvature0, deltaL, EI);
dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
curvature0, deltaL, EI);
f(1:6) = f(1:6) + dF;
J(1:6,1:6) = J(1:6,1:6) + dJ;

% Viscous force
f = f + C * ( q - q0 ) / dt;
J = J + C / dt;

% Weight
f = f - W;% No need to update the Jacobian

% Update
q = q - J \ f;% Main part of Newton's method
err = sum( abs(f) );
end
% Update
u = (q - q0) / dt; % Velocity
q0 = q; % Old position
%figure(1);
figure(counter+1);
if c==node_number(counter)
xcoords=q(1:2:end) %x1,x2,x3
ycoords=q(2:2:end) %y1,y2,y3
plot( q(1:2:end), q(2:2:end), 'ro-');
if(node_number(counter)~=999)
counter=counter+1;
end
axis equal;
xlabel('Position, m [meter]');
ylabel('Position, m [meter]');
end
% Store
all_mid_y(c) = q(4);% appending the position of the middle node as the simulation is running
all_mid_v(c) = u(4);% appending the velocity of the middle node as the simulation is running
end
figure(counter+2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');
figure(counter+3);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_y, 'k-');
xlabel('Time, t [sec]');
ylabel('Position of mid-node, v [meter/sec]');