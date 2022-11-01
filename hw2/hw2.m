%%
clc;
clear all;
close all;

%% Global variables
global m dt unconsInd tol ScaleSolver maximum_iter
global Fg
global EI GJ voronoiRefLen kappaBar EA refLen
%% Parameters inputs 
vertices = 100; % number of vertices
dt = 0.01; % Time step
Rod_L = 0.2; % rod length
rho = 1000; % density
r0 = 1e-3; % cross sectional radius of the rod
Y = 10e6; % Young's modulus
nu = 0.5; % Poisson ratio
G = Y / (2*(1+nu)); % Shear modulus
natR = 0.02; % Natural radius of curvature
g = [0;0;-9.81];% gravity
tol = 1e-3; % Tolerance on force function (to be ...
% multiplied by characteristic bending force)
maximum_iter = 100; % maximum iterations in Newton's solver
totalTime = 5; % total simulation time
% Stiffness properties
EI = Y * pi * r0^4 / 4; % bending stiffness
EA = Y * pi * r0^2; % stretching stiffness
GJ = G * pi * r0^4 / 2; % twisting stiffness

num_dof = 4*vertices - 1; % number of degrees of freedom
edges =  vertices - 1; % number of edges
dm = (pi * r0^2 * Rod_L) * rho / edges;
ScaleSolver = EI / Rod_L^2; % Characteristic bending force
%% Mass
m = zeros(num_dof, 1);
for c=1:vertices
    if c==1 || c==vertices
        m( 4*(c-1) + 1: 4*(c-1) + 3 ) = dm/2;
    else
        m( 4*(c-1) + 1: 4*(c-1) + 3 ) = dm;        
    end        
end
for c=1:edges
    m( 4*c ) = dm/2 * r0^2; % I = 1/2*m*r^2
end

%% Gravity
garr = zeros(num_dof, 1);
for c=1:vertices
    garr( 4*(c-1) + 1: 4*(c-1) + 3 ) = g;
end
Fg = m .* garr; % Weight vector
%% Geometry of the rod
nodes = zeros(vertices, 3);
dTheta = (Rod_L/edges) / natR;
for c=1:vertices
    nodes(c,1) = natR * cos( (c-1) * dTheta );
    nodes(c,2) = natR * sin( (c-1) * dTheta );
end
%% Reference length
refLen = zeros(edges, 1);
for c=1:edges
    dx = nodes(c+1,:) - nodes(c,:);
    refLen(c) = norm(dx);
end
%% Voronoi length
voronoiRefLen = zeros(vertices, 1);
for c=1:vertices
    if c==1
        voronoiRefLen(c) = 0.5*refLen(c);
    elseif c==vertices
        voronoiRefLen(c) = 0.5*refLen(c-1);
    else
        voronoiRefLen(c) = 0.5*(refLen(c-1) + ...
            refLen(c));
    end
end

%% Reference Frames 
A1 = zeros(edges, 3); % First Reference director, u
A2 = zeros(edges, 3); % Second Reference director, v

tangent = zeros(edges, 3); % Tangent on the first edge 
for c=1:edges
    dx = nodes(c+1,:) - nodes(c,:);
    tangent(c,:) = dx / norm(dx);
end

% Figure out a good choice for d1(1,:)
t0 = tangent(1,:);
t1 = [0 0 -1];
d1Tmp = cross(t0, t1);
if abs(d1Tmp) < 1e-6
    t1 = [ 0 1 0 ];
    d1Tmp = cross(t0, t1);
end
A1(1,:) = d1Tmp / norm(d1Tmp);
d2Tmp = cross(t0, A1(1,:));
A2(1,:) = d2Tmp / norm(d2Tmp);

for c=2:edges
    t0 = tangent(c-1,:); %tangent on the c-1 th edge 
    t1 = tangent(c,:); % tangent on the cth edge 
    A1_old = A1(c-1,:);
    A1_new = parallel_transport(A1_old, t0, t1);
    A1_new = A1_new / norm(A1_new);
    A1(c,:) = A1_new;
    A2_new = cross(t1, A1_new);
    A2(c,:) = A2_new / norm(A2_new);
end

%% Old and new degrees of freedom vector
x0 = zeros(num_dof, 1); % Old dof (i.e. q0)
for c=1:vertices
    x0( 4*(c-1) + 1 ) = nodes(c,1); % x-coord
    x0( 4*(c-1) + 2 ) = nodes(c,2); % y-cood
    x0( 4*(c-1) + 3 ) = nodes(c,3); % z-coord
end
x0(4:4:end) = 0; % theta

x = x0; % New dof (i.e. q)
u = (x - x0) / dt; % Old velocity vector

%% Constrained (fixed) and unconstrained (free) indices
consInd = 1:7; % Clamped; fixed_index
unconsInd = 8:num_dof; % free_index

%% Material director
theta = x(4:4:end);
[m1, m2] = computeMaterialDirectors(A1, A2, theta);

%% Reference twist
refTwist = zeros(vertices, 1); % Reference twist 
refTwist = getRefTwist( A1, tangent, refTwist );

%% Natural curvature
kappaBar = getkappa( x, m1, m2 );


%% Simulation loop
Nsteps = round(totalTime/dt); % number of time steps

endZ = zeros(Nsteps, 1);

ctime = 0; % Current time

A1_old = A1;
A2_old = A2;
refTwist_old = refTwist;

for timeStep = 1:Nsteps % Main simulation step
    
    fprintf('t = %f\n', ctime);
    
    [x, A1, A2, refTwist] = objfun(x0, u, A1_old, A2_old, refTwist_old); % Main function that steps forward one time step

    % Plot
    x_coord = x(1:4:end);
    y_coord = x(2:4:end);
    z_coord = x(3:4:end);
    figure(1);
    clf();
    plot3(x_coord, y_coord, z_coord, 'ko-');
    title('Simulation at Time=5 secs');
    axis equal
    xlabel('x [meters]');
    ylabel('y [meters]');
    zlabel('z [meters]');
    endZ(timeStep) = x(end);

    u = (x - x0) / dt; % Update velocity
    x0 = x; % New position becomes old position
    A1_old = A1; % New reference director becomes old reference director
    A2_old = A2;
    refTwist_old = refTwist; % New reference twist becomes old reference twist
    
    
    ctime = ctime + dt; % Current time
end

%% plotting the results 
figure(2);
time_array = (1:Nsteps) * dt;
plot( time_array, endZ, 'ko-');
title('End Z coordinate vs Time for N = 100');
xlabel('Time, t [seconds]');
ylabel('z-coord of last node, \delta_z [meters]');