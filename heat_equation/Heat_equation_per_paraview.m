clear all
clc
close all

% Generate Geometry
%------------------------------------
display_nurbs_surface_square;
pause( 0.5 );
    
nurbs_initial=nurbs;

% Choice of the NURBS degree, regularity of the basis and mesh size
p_vector = [2]; % non so quali siano adatti al nostro problema, ho messo valori un po' a caso
n_elem = 2.^[ 5 ]; % in one parametric direction
h_vector = 1 ./ n_elem; % characteristic mesh size for the geometry

    
% set phyiscal DATA, boundary conditions
%----------------------------------------
mu = @( x, y )( 1e5 + 0 * x );   % diffusion coefficient
b = @( x, y ) ( [  0 * x.*y; 0 * x.*y ] );
f = @( x, y )( (8*pi^2-1).*sin(2*pi*x).*sin(2*pi*y));    % source term
f_time = @( t ) (exp(-t) );  
rho = @( x, y )( 1 + 0* x.*y );   % coefficient for time dependent term
    
drchlt_sides = [ 1 2 3 4 ];   % indexes for Dirichlet faces
nmnn_sides   = [];   % indexes for Neuman faces
    
g = @( x, y, ind )( sin(2*pi*x).*sin(2*pi*y).*exp(-ind)  ); % Dirichlet BCs
h = @( x, y, ind )( 0 + 0 * x );    % Neuman data (for different boundaries)
    
drchlt_imposition_type = 'int'; % 'L2' = L2 projection, 
                                % 'int' = interpolation at control points

% set exact solution, if available
%-----------------------------------
uex = @(x,y,t) sin(2*pi*x).*sin(2*pi*y).*exp(-t);
graduex = @(x,y,t) [2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-t); 2*pi*cos(2*pi*y).*sin(2*pi*x).*exp(-t)];

% time and time step
%--------------------
Tf = 0.2;
dt = 0.001;
Nt = round( ( Tf - 0 ) / dt );
                                    
% parameter for theta method
%----------------------------
theta = 0.5;
                   
% output settings
%------------------
output_file_name = 'TestT2D_';
n_pts_viz = 65;  % number of points for visualization
vtk_pts = { linspace( knots{ 1 }( 1 ), knots{ 1 }( end ), n_pts_viz ), ...
            linspace( knots{ 2 }( 1 ), knots{ 2 }( end ), n_pts_viz ) };   
                                                    
                               
figure( 25 )
nrbkntplot( nurbs );
view( 0, 90 );
pause( 0.5 );

p=p_vector(1);
k=1;
nurbs = nrbdegelev( nurbs_initial, [ ( p - 2 ) ( p - 2 ) ]);
           
new_knots_notrepeated = [ 1 : ( n_elem - 1 ) ] / n_elem;
new_knots = [];

for m = p - 1 : -1 : k
    new_knots = [ new_knots, new_knots_notrepeated ];
end

new_knots = sort( new_knots );
nurbs = nrbkntins( nurbs, { new_knots, new_knots } );

ndof_tot_vector = [nurbs.number( 1 ) * nurbs.number( 2 ) ];

           
%==========================================================================

% Generate Space & mesh infos
%-----------------------------
geometry = geo_load( nurbs );
knots = geometry.nurbs.knots;
[ qn, qw ] = msh_set_quad_nodes( knots, ...
                                 msh_gauss_nodes( geometry.nurbs.order ) );
msh = msh_2d( knots, qn, qw, geometry );
space = sp_nurbs_2d( geometry.nurbs, msh );

% assemble the matrix and vectors whose data do not depend on time.
%-----------------------------------------------------------------
matrix_A = op_gradu_gradv_tp( space, space, msh, mu );
matrix_M = op_u_v_tp( space, space, msh, rho );

rhs_F = op_f_v_tp( space, msh, f );
rhs_Neuman = zeros( size( rhs_F ) );

% Set Neuman BCs
%------------------
for iside = nmnn_sides
  msh_side = msh_eval_boundary_side( msh, iside );
  sp_side  = sp_eval_boundary_side( space, msh_side );

  x = squeeze( msh_side.geo_map( 1, :, : ) );
  y = squeeze( msh_side.geo_map( 2, :, : ) );
  hval = reshape( h( x, y, iside ), msh_side.nqn, msh_side.nel );

  rhs_Neuman( sp_side.dofs ) = rhs_Neuman( sp_side.dofs ) ...
                        + op_f_v( sp_side, msh_side, hval );
end

% Set Dirichlet BCs
%----------------------
switch drchlt_imposition_type
    
    case 'L2'
        % L2 projection
        [ u_drchlt, drchlt_dofs ] = sp_drchlt_l2_proj( space, msh, ...
                                               g, drchlt_sides );
    case 'int'
        % Interpolation at control points
        drchlt_dofs_iside = [];
        u_drchlt_iside = [];
        for iside = drchlt_sides
            drchlt_dofs_iside = [ drchlt_dofs_iside, space.boundary( iside ).dofs ];
            switch iside
                case 1
                    i1_iside = 1;
                    i2_iside = 1 : nurbs.number( 2 );
                    n_dof_iside = nurbs.number( 2 );
                case 2
                    i1_iside = nurbs.number( 1 );
                    i2_iside = 1 : nurbs.number( 2 );
                    n_dof_iside = nurbs.number( 2 );
                case 3
                    i1_iside = 1 : nurbs.number( 1 );
                    i2_iside = 1;
                    n_dof_iside = nurbs.number( 1 );
                case 4
                    i1_iside = 1 : nurbs.number( 1 );
                    i2_iside = nurbs.number( 2 );
                    n_dof_iside = nurbs.number( 1 );
            end
            x_iside = reshape( nurbs.coefs( 1, i1_iside, i2_iside ) ./ nurbs.coefs( 4, i1_iside, i2_iside ), n_dof_iside, 1 );
            y_iside = reshape( nurbs.coefs( 2, i1_iside, i2_iside ) ./ nurbs.coefs( 4, i1_iside, i2_iside ), n_dof_iside, 1 );    
            u_drchlt_iside = [ u_drchlt_iside; g( x_iside, y_iside, iside ) ];       
        end
        [ drchlt_dofs, i_drchlt_dofs ] = unique( drchlt_dofs_iside ); 
        u_drchlt = u_drchlt_iside( i_drchlt_dofs );
        
end
    %
    int_dofs = setdiff( 1 : space.ndof, drchlt_dofs );
    
    
    % set intial solution
    %---------------------
    % Definisci la funzione di inizializzazione
    u_init = @(x,y) sin(2*pi*x).*sin(2*pi*y);
    
    coefs = reshape(nurbs.coefs, [], nurbs.number(1) * nurbs.number(2));
    X = coefs(1,:) ./ coefs(4,:);
    Y = coefs(2,:) ./ coefs(4,:);
    % Valuta la funzione iniziale su queste coordinate
    u_init_values = arrayfun(@(x, y) sin(2*pi*x).*sin(2*pi*y), X, Y);
    % Assegna questi valori iniziali a u
    u = zeros(space.ndof, 1);
    u = u_init_values';

    % Loop over time
    %==========================================
    for n = 0 : Nt - 1
        
        u_old = u;
        
        time_n = n * dt;
        time_np1 = ( n + 1 ) * dt;
        
        rhs_time_n = f_time( time_n ) * rhs_F + rhs_Neuman;
        rhs_time_np1 = f_time( time_np1 ) * rhs_F + rhs_Neuman;
        
        matrix = matrix_M + theta * dt * matrix_A;
        rhs = ( matrix_M - ( 1 - theta ) * dt * matrix_A ) * u_old ...
              + dt * ( theta * rhs_time_np1 + ( 1 - theta ) * rhs_time_n );
          
        % apply Dirichlet BCs
        %---------------------
        u = zeros( space.ndof, 1 );
        u( drchlt_dofs ) = u_drchlt;
        rhs( int_dofs ) = rhs(int_dofs) ...
                      - matrix( int_dofs, drchlt_dofs ) * u_drchlt;
    
    
        % Solve the linear system (direct method)
        %------------------------------------------
        u( int_dofs ) = matrix( int_dofs, int_dofs ) \ rhs( int_dofs );
        

        % save vtk file
        %---------------%
        % SALVATAGGIO DEI DATI IN PIU' FILE
        % NOME DELLA CARTELLA DA CAMBIARE OGNI VOLTA
        output_folder = 'results/results_a';
        
        % Crea la cartella se non esiste
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        
        % NOME DEL FILE DA CAMBIARE OGNI VOLTA
        output_file_name_n = sprintf('%s/results_a_%04d.vtk', output_folder, n);
        sp_to_vtk(u(:, end), space, geometry, vtk_pts, output_file_name_n, 'u');
        
        fprintf('The result is saved in the file: %s \n', output_file_name_n);
        
           
    end

