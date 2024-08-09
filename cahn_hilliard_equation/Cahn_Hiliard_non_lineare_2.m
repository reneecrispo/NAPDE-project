clear all;
clc;
close all;

% Generate Geometry
%------------------------------------
%
display_nurbs_surface_square;
pause(0.5);

nurbs_initial = nurbs;

% Choice of the NURBS degree, regularity of the basis and mesh size
p_vector = [2]; 
n_elem = 2^4; % in one parametric direction
h = 1 / n_elem; % characteristic mesh size for the geometry

% set physical DATA, boundary conditions
%----------------------------------------
mu = @(x,y) 1;   % bilaplacian coefficient
rho = @(x, y) (1 + 0 * x .* y);   % coefficient for time dependent term

drchlt_sides = [1 2 3 4];   % indexes for Dirichlet faces
nmnn_sides = [];   % indexes for Neuman faces


drchlt_imposition_type = 'int'; % 'L2' = L2 projection, 
                                % 'int' = interpolation at control points

% time and time step
%--------------------
Tf = 0.1;
dt = 0.001;
Nt = round(Tf / dt);

% parameter for theta method
%----------------------------
theta = 0;

% output settings
%------------------
output_file_name = 'TestT2D_';
n_pts_viz = 65;  % number of points for visualization
knots = nurbs_initial.knots; % aggiunto
vtk_pts = {linspace(knots{1}(1), knots{1}(end), n_pts_viz), ...
           linspace(knots{2}(1), knots{2}(end), n_pts_viz)};

% parte che crea errori                               
figure(25);
nrbkntplot(nurbs); 
view(0, 90);
pause(0.5);

p = p_vector(1); % definisci p
k = 1; % definisci k
nurbs = nrbdegelev(nurbs_initial, [(p - 2) (p - 2)]);

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

% matrix_T è il tensore associato all'integrale del prodotto di tre
% funzioni di base 

% Numero di funzioni di base
N = size(matrix_M, 1);

% Inizializzare la matrice cubica
matrix_T = zeros(N, N, N);

% Calcoliamo la matrice cubica usando il prodotto di Kronecker
for i = 1:N
    for j = 1:N
        % prodotto di Kronecker delle righe
        kron_prod = kron(matrix_M(i, :), matrix_M(j, :));
       
        % Riorganizziamo il risultato in una matrice [N, N]
        kron_matrix = reshape(kron_prod, [N, N]);
        
        % Salviamo la matrice risultante nella matrice cubica
        matrix_T(i, j, :) = kron_matrix * matrix_M(j, :)';
    end
end

% Set Dirichlet BCs
%----------------------
switch drchlt_imposition_type
    case 'L2'
        % L2 projection
        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj(space, msh, g, drchlt_sides);
    case 'int'
        % Interpolation at control points
        drchlt_dofs_iside = [];
        u_drchlt_iside = [];
        for iside = drchlt_sides
            drchlt_dofs_iside = [drchlt_dofs_iside, space.boundary(iside).dofs];
            switch iside
                case 1
                    i1_iside = 1;
                    i2_iside = 1:nurbs.number(2);
                    n_dof_iside = nurbs.number(2);
                case 2
                    i1_iside = nurbs.number(1);
                    i2_iside = 1:nurbs.number(2);
                    n_dof_iside = nurbs.number(2);
                case 3
                    i1_iside = 1:nurbs.number(1);
                    i2_iside = 1;
                    n_dof_iside = nurbs.number(1);
                case 4
                    i1_iside = 1:nurbs.number(1);
                    i2_iside = nurbs.number(2);
                    n_dof_iside = nurbs.number(1);
            end
            x_iside = reshape(nurbs.coefs(1, i1_iside, i2_iside) ./ nurbs.coefs(4, i1_iside, i2_iside), n_dof_iside, 1);
            y_iside = reshape(nurbs.coefs(2, i1_iside, i2_iside) ./ nurbs.coefs(4, i1_iside, i2_iside), n_dof_iside, 1);    
            % u_drchlt_iside = [u_drchlt_iside; g(x_iside, y_iside, iside)];       
        end
        [drchlt_dofs, i_drchlt_dofs] = unique(drchlt_dofs_iside); 
        % u_drchlt = u_drchlt_iside(i_drchlt_dofs);
end

int_dofs = setdiff(1:space.ndof, drchlt_dofs);

% initial condition
%=========================================
u = rand(space.ndof, 1);

% Loop over time
%==========================================
for n = 0:Nt-1
    u_old = u;

    time_n = n * dt;
    time_np1 = (n + 1) * dt;

    % passaggi per calcolare u_old' * matrix_T * u_old
    B = squeeze(sum(matrix_T .* reshape(u_old, [1, 1, numel(u_old)]), 3));
    result1 = B .* u_old';
    result2 = B .* (u_old .* u_old)';

    u = matrix_M \ (matrix_M * u_old - dt * matrix_A * u_old - dt * result2.^2 * (u_old .* u_old) + dt * result1 * u_old);

end
