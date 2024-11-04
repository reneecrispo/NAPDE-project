clear all;
clc;
close all;


% Generate Geometry
%------------------------------------
display_nurbs_surface_square;
%pause(0.5);

nurbs_initial = nurbs;

p_vector = [2]; 
n_elem = 2^4; 
h = 1 / n_elem; 

% set physical DATA, boundary conditions
%----------------------------------------
mu = @(x,y) 0.4*h^2;   % bilaplacian coefficient
rho = @(x,y) (1 + 0 * x .* y);   % coefficient for time dependent term


% time and time step
%--------------------
dt = 1e-3;
Tf = 800*dt;
Nt = round(Tf / dt);


% output settings
%------------------
output_file_name = 'TestT2D_';
n_pts_viz = 129;  
knots = nurbs_initial.knots;
vtk_pts = {linspace(knots{1}(1), knots{1}(end), n_pts_viz), ...
           linspace(knots{2}(1), knots{2}(end), n_pts_viz)};

p = p_vector(1); 
k = 1;
nurbs = nrbdegelev(nurbs_initial, [(p - 1) (p - 1)]);

new_knots_notrepeated = [1:(n_elem - 1)] / n_elem;
new_knots = [];

for m = p - 1:-1:k
    new_knots = [new_knots, new_knots_notrepeated];
end

new_knots = sort(new_knots);
nurbs = nrbkntins(nurbs, {new_knots, new_knots});

ndof_tot_vector = [nurbs.number(1) * nurbs.number(2)];

% Generate Space & mesh infos
%-----------------------------
geometry = geo_load(nurbs);
knots = geometry.nurbs.knots;
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(geometry.nurbs.order));
msh = msh_2d(knots, qn, qw, geometry);
space = sp_nurbs_2d(geometry.nurbs, msh);

% assemble the matrix and vectors whose data do not depend on time.
%-----------------------------------------------------------------
matrix_A = op_laplaceu_laplacev_tp(space, space, msh, mu);
matrix_M = op_u_v_tp(space, space, msh, rho);

% Dirichlet homogeneous boudary conditions on a frame of width of 2 DOF
%---------------------------------------------------------------------------------------
drchlt_dofs = []; 

N = sqrt(space.ndof);
u = zeros(space.ndof, 1);
u_drchlt = 0; 

% Definitions of DOFs
drchlt_dofs = unique([1:N, ...                    % First row
                      N*(N-1)+1:N^2, ...          % Last row
                      1:N:N*(N-1)+1, ...          % First column
                      N:N:N^2, ...                % Last column
                      N+1:2*N, ...                % Second row
                      N*(N-2)+1:N*(N-1), ...      % Second-last row
                      2:N:N^2-N+2, ...            % Second column
                      N-1:N:N^2-1]);              % Second-last column

u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff(1:space.ndof, drchlt_dofs);


% initial condition
%=========================================
coefs = reshape(nurbs.coefs, [], nurbs.number(1) * nurbs.number(2));
X = coefs(1,:) ;
Y = coefs(2,:) ; 


for i = 1:length(X)
    if((X(1,i) >= 0.50 && X(1,i) <= 0.60) && ((Y(1,i) >= 0.50 && Y(1,i) <= 0.55)))
     u_init_values(i)= 1.0;
        else 
       u_init_values(i) = 0.0;
    end
end

u = zeros(space.ndof, 1);
u(1:length(u_init_values)) = u_init_values;

% Nutrient definition:
% --------------------------------------------------------------
nut_g = zeros(space.ndof, 1); 

% Initialize an array for the boundary degrees of freedom
drchlt_dofs_nut = [];

% Loop through the sides of the boundary
for iside = 1:numel(space.boundary)
    drchlt_dofs_nut = union(drchlt_dofs_nut, space.boundary(iside).dofs);  % Union the dofs of each side
end

% Set the values to 1 at the boundary nodes
nut_g(drchlt_dofs_nut) = 1;

Dn = @(x,y) (1 + 0 * x .* y); 
matrix_A_nut = op_gradu_gradv_tp(space, space, msh, Dn);
Fg = -matrix_A_nut * nut_g;

% Set the rows and columns of matrix_A_nut associated with Dirichlet degrees of freedom to zero
for i = 1:length(drchlt_dofs_nut)
    dof = drchlt_dofs_nut(i);
    
    % Set the corresponding row and column for dof in matrix_A_nut to zero
    matrix_A_nut(dof, :) = 0;
    matrix_A_nut(:, dof) = 0;
    
    % Set the diagonal element to 1 to maintain non-singularity
    matrix_A_nut(dof, dof) = 1;
    
    % Set the corresponding term in Fg to zero
    Fg(dof) = 0;
end

% matrix_M_nut = op_u_v_tp(space, space, msh, rho);
nut_dot = zeros(space.ndof, 1);
nut = zeros(space.ndof, 1);
nut_dot = matrix_A_nut \ Fg;
nut = nut_dot + nut_g;

% nut(drchlt_dofs_nut) = 1;
P0 = 5;


% Loop over time
%==========================================
for n = 0 : (Nt-1)
    u_old = u(:, end);

    time_n = n * dt;
    time_np1 = (n + 1) * dt;

    g = @(u) 2*(- 2 - 4*u.^2 + 6*u);    
    dg = @(u) 2*(-8*u + 6);       
     
    matrix_T = @(x) op_u_v_tp_cahn_hilliard_non_lin(space, space, msh,g,x);
    matrix_derT = @(x) op_u_v_tp_cahn_hilliard_non_lin(space, space, msh,dg,x);

    fun = @(x) (matrix_M + dt * matrix_A - dt * P0*matrix_M.*nut - dt * matrix_T(x)) * x ...
            - matrix_M * u_old;
    
    J = @(x) matrix_M + dt * matrix_A - dt *P0* matrix_M.*nut - dt * matrix_T(x) ...
           - dt * matrix_derT(x) * x * ones(size(x))' * matrix_M';
   
    niter = 10;
    toll = 1e-3;

    [u, it] = newtonsys(u_old, niter, toll, fun, J);
    u(drchlt_dofs) = 0; 
    it;

    output_folder = 'results/results_CN_h';
    output_file_name_n = sprintf('%s/results_CN_h_%04d', output_folder, n+1);

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    sp_to_vtk(u(:, end), space, geometry, vtk_pts, output_file_name_n, 'u');

    fprintf('The result is saved in the file: %s \n', output_file_name_n);
end