clear all;
clc;
close all;

%format long 

% Generate Geometry
%------------------------------------
%
display_nurbs_surface_square;
%pause(0.5);

nurbs_initial = nurbs;

% Choice of the NURBS degree, regularity of the basis and mesh size
p_vector = [2]; 
n_elem = 2^4; % in one parametric direction
h = 1 / n_elem; % characteristic mesh size for the geometry

% set physical DATA, boundary conditions
%----------------------------------------
mu = @(x,y) 0.3*h^2;   % bilaplacian coefficient
rho = @(x,y) (1 + 0 * x .* y);   % coefficient for time dependent term

drchlt_sides = [];   % indexes for Dirichlet faces
nmnn_sides = [1 2 3 4];   % indexes for Neumann faces

drchlt_imposition_type = 'int'; % 'L2' = L2 projection, 
                                % 'int' = interpolation at control points

% time and time step
%--------------------
dt = 1e-4;
Tf = 1000*dt;
Nt = round(Tf / dt);

% parameter for theta method
%----------------------------
theta = 0;

% output settings
%------------------
output_file_name = 'TestT2D_';
n_pts_viz = 129;  % number of points for visualization
knots = nurbs_initial.knots; % aggiunto
vtk_pts = {linspace(knots{1}(1), knots{1}(end), n_pts_viz), ...
           linspace(knots{2}(1), knots{2}(end), n_pts_viz)};


p = p_vector(1); % definisci p
k = 1; % definisci k
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

% Impostazione delle condizioni di Dirichlet omogenee su una cornice di spessore due DOF
%---------------------------------------------------------------------------------------
drchlt_dofs = []; % Lista dei DOFs per la cornice

% Impostazione della cornice con condizioni di Dirichlet su una griglia 19x19
N = sqrt(space.ndof); % Numero di nodi per lato (19x19 griglia)
u = zeros(space.ndof, 1); % Vettore dei valori su tutta la griglia
u_drchlt = 0; % Condizione Dirichlet omogenea (zero)

% Definizione dei DOFs per la cornice doppia (primo e secondo strato di nodi ai bordi)
drchlt_dofs = unique([1:N, ...                    % Prima riga
                      N*(N-1)+1:N^2, ...          % Ultima riga
                      1:N:N*(N-1)+1, ...          % Prima colonna
                      N:N:N^2, ...                % Ultima colonna
                      N+1:2*N, ...                % Seconda riga
                      N*(N-2)+1:N*(N-1), ...      % Penultima riga
                      2:N:N^2-N+2, ...            % Seconda colonna
                      N-1:N:N^2-1]);            % Penultima colonna
% Imposta il valore di Dirichlet (zero) sui DOFs della cornice
u(drchlt_dofs) = u_drchlt;

% Visualizza il numero di DOFs e la conferma dell'imposizione Dirichlet
fprintf('Numero di DOFs nella cornice doppia: %d\n', length(drchlt_dofs));
disp('Condizioni di Dirichlet imposte con successo sui DOFs della cornice doppia.');
% Trova i gradi di libertà interni
int_dofs = setdiff(1:space.ndof, drchlt_dofs);


% initial condition
%=========================================
coefs = reshape(nurbs.coefs, [], nurbs.number(1) * nurbs.number(2));
X = coefs(1,:) ; % Coordinate X dei nodi
Y = coefs(2,:) ; % Coordinate Y dei nodi


for i = 1:length(X)
    if((X(1,i) >= 0.30 && X(1,i) <= 0.70) && ((Y(1,i) >= 0.45 && Y(1,i) <= 0.55)))
      %u_init_values(i) = 0.5 + 0.2 * sin(1*pi*X(1,i)).^0.1 .* sin(1*pi*Y(1,i)).^0.1; 
       % u_init_values(i)= 0.2 * cos(8*pi*X(1,i)).^0.1 .* cos(8*pi*Y(1,i)).^0.1;
       u_init_values(i)= 1;
        %u_init_values(i)= 0.2+ 0.1 * cos(18*pi*X(1,i)).* cos(18*pi*Y(1,i));
    else 
       u_init_values(i) = 0;
    end
end


% Assegno i valori iniziali di u allo spazio dei gradi di libertà
u = zeros(space.ndof, 1);
u(1:length(u_init_values)) = u_init_values;

% Definizione dei nutrienti:
% --------------------------------------------------------------
P_0 = 10;
gamma_n = @(x,y) P_0 * (1 - 0.5*x + 0*y);

matrix_M_tilde = op_u_v_tp(space, space, msh, gamma_n);


% Loop over time
%==========================================
for n = 0 : (Nt-1)
    u_old = u(:, end);

    time_n = n * dt;
    time_np1 = (n + 1) * dt;

    % g è la funzione g(u) con raccolta una u 
    g = @(u) 2*(- 2 - 4*u.^2 + 6*u);    
    dg = @(u) 2*(-8*u + 6);
    %GG = @(u) 2*(- 2 - 4*u.^2 + 6*u) .* u;
    %dGG = @(u) 2*(- 2 - 4*u.^2 + 6*u) + 2*(-8*u + 6) .* u;
       
     
    matrix_T = @(x) op_u_v_tp_cahn_hilliard_non_lin(space, space, msh,g,x);
    matrix_derT = @(x) op_u_v_tp_cahn_hilliard_non_lin(space, space, msh,dg,x);
    %matrix_derT2 = @(x) op_u_v_tp_cahn_hilliard_non_lin(space, space, msh,dGG,x);

    fun = @(x) (matrix_M + dt * matrix_A - dt * matrix_M_tilde - dt * matrix_T(x)) * x ...
            - matrix_M * u_old;
   
    J = @(x) matrix_M + dt * matrix_A - dt * matrix_T(x) ...
       - dt * matrix_derT(x) * x * ones(size(x))' * matrix_M';
   
    niter = 15;
    toll = 1e-3;

    [u, it] = newtonsys(u_old, niter, toll, fun, J);
    u(drchlt_dofs) = 0; %imposizione condizioni di Dirichelet
    it;

    % Save results to files
    output_folder = 'results/results_CP_o';
    % File name for each time step
    output_file_name_n = sprintf('%s/results_CP_o_%04d', output_folder, n+1);

    % Create the folder if it does not exist
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % File name for each time step
    %output_file_name_n = sprintf('%s/results_ee_%04d.vtk', output_folder, n);

    sp_to_vtk(u(:, end), space, geometry, vtk_pts, output_file_name_n, 'u');

    fprintf('The result is saved in the file: %s \n', output_file_name_n);
end