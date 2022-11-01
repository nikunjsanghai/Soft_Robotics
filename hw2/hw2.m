%%
clc;
clear all;
close all;

%% Global variables
global m time_step unconsInd tol SSol max_iter
global Fg
global EI GJ v_RefLen kappaBar EA refLen
%% Parameters
vertices = 100; % number of vertices
time_step = 0.01; % Time step
Rod_Len = 0.2; % rod length
rho = 1000; % density
r0 = 1e-3; % cross sectional radius of the rod
Y = 10e6; % Young's modulus
nu = 0.5; % Poisson ratio
G = Y / (2*(1+nu)); % Shear modulus
nat_Radius = 0.02; % Natural radius of curvature
g = [0;0;-9.81];
tol = 1e-3; % Tolerance on force function 
max_iter = 100; % maximum iterations in Newton's solver
total_T = 5; % total simulation time

% Stiffness properties
EI = Y * pi * r0^4 / 4; % bending stiffness
EA = Y * pi * r0^2; % stretching stiffness
GJ = G * pi * r0^4 / 2; % twisting stiffness

% few other parameters
ndof = 4*vertices - 1; % number of degrees of freedom
edges =  vertices - 1; % number of edges
dm = (pi * r0^2 * Rod_Len) * rho / edges;
SSol = EI / Rod_Len^2; % Characteristic bending force
%% Geometry of the rod
nodes = zeros(vertices, 3);
dTheta = (Rod_Len/edges) / nat_Radius;
for c=1:vertices
    nodes(c,1) = nat_Radius * cos( (c-1) * dTheta );
    nodes(c,2) = nat_Radius * sin( (c-1) * dTheta );
end
%% Mass
m = zeros(ndof, 1);
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
garr = zeros(ndof, 1);
for c=1:vertices
    garr( 4*(c-1) + 1: 4*(c-1) + 3 ) = g;
end
Fg = m .* garr; % Weight vector
%% Reference Length(Edge Length)
refLen = zeros(edges, 1);
for c=1:edges %loop over the edges 
    dx = nodes(c+1,:) - nodes(c,:);
    refLen(c) = norm(dx);
end

%% Voronoi length (Length associated with each node)
v_RefLen = zeros(vertices, 1);
for c=1:vertices %loop over nodes
    if c==1
        v_RefLen(c) = 0.5*refLen(c);
    elseif c==vertices
        v_RefLen(c) = 0.5*refLen(c-1);
    else
        v_RefLen(c) = 0.5*(refLen(c-1) + ...
            refLen(c));
    end
end


%% Reference Frame
d1 = zeros(edges, 3); % Reference director, u
d2 = zeros(edges, 3); % Reference director, v

tangent = zeros(edges, 3);
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
d1(1,:) = d1Tmp / norm(d1Tmp);
d2Tmp = cross(t0, d1(1,:));
d2(1,:) = d2Tmp / norm(d2Tmp);

for c=2:edges
    t0 = tangent(c-1,:);
    t1 = tangent(c,:);
    d1_old = d1(c-1,:);
    A1_new = parallel_transport(d1_old, t0, t1);
    A1_new = A1_new / norm(A1_new);
    d1(c,:) = A1_new;
    A2_new = cross(t1, A1_new);
    d2(c,:) = A2_new / norm(A2_new);
end

%% Old and new degrees of freedom vector
x0 = zeros(ndof, 1); % Old dof (i.e. q0)
for c=1:vertices
    x0( 4*(c-1) + 1 ) = nodes(c,1); % x-coord
    x0( 4*(c-1) + 2 ) = nodes(c,2); % y-cood
    x0( 4*(c-1) + 3 ) = nodes(c,3); % z-coord
end
x0(4:4:end) = 0; % theta

x = x0; % New dof (i.e. q)
u = (x - x0) / time_step; % Old velocity vector

%% Constrained (fixed) and unconstrained (free) indices
consInd = 1:7; % Clamped; fixed_index
unconsInd = 8:ndof; % free_index

%% Material Frame
theta = x(4:4:end);
[m1, m2] = computeMaterialDirectors(d1, d2, theta);

%% Reference twist
refTwist = zeros(vertices, 1);
refTwist = getRefTwist( d1, tangent, refTwist );

%% Natural curvature
kappaBar = getkappa( x, m1, m2 );


%% Simulation 
Nsteps = round(total_T/time_step); % number of time steps

endZ = zeros(Nsteps, 1);

ctime = 0; % Current time

d1_old = d1;
d2_old = d2;
refTwist_old = refTwist;

for timeStep = 1:Nsteps % Main simulation step
    
    fprintf('t = %f\n', ctime);
    
    [x, d1, d2, refTwist] = objfun(x0, u, d1_old, d2_old, refTwist_old); % function that Iterates one time step

    % Plot
    x_coord = x(1:4:end);
    y_coord = x(2:4:end);
    z_coord = x(3:4:end);
    figure(1);
    clf();
    plot3(x_coord, y_coord, z_coord, 'co-');
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    endZ(timeStep) = x(end);

    u = (x - x0) / time_step; % Update velocity
    x0 = x; % New position becomes old position
    d1_old = d1; % New reference director becomes old reference director
    d2_old = d2;
    refTwist_old = refTwist; % New reference twist becomes old reference twist
    
    
    ctime = ctime + time_step; % Current time
end

%%
figure(2);
tarray = (1:Nsteps) * time_step;
plot( tarray, endZ, 'go-');
title('End Z coordinate vs Time for N = 100');
xlabel('Time, t [s]');
ylabel('z-coord of last node, \delta_z [m]');