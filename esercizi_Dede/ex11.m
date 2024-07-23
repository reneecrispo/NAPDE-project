%--------------------------------------------------------------------------
% Solve the time dependent advection-diffusion equation in 2D by 
% NURBS-based IGA, Galerkin method and the theta-method. 
% Based on GeoPDEs. Imposition of the Dirichlet boundary conditions 
% by means of L2 projection or "interpolation at the control points".
%
%--------------------------------------------------------------------------
%
% \rho du/dt - \nabla \cdot (\mu \nabla u) + b \cdot \nabla u = f,  in \Omega x (0,T)
%                               u = g,    on \Gamma_Dirichlet x (0,T)
%            \mu \nabla u \cdot n = h,    on \Gamma_Neumann x (0,T)
%                           u( t = 0 ) = u0,   in \Omega
% 
% f = f(t,x)
% \mu = \mu(x), \rho=\rho(x), b = b(x), g=g(x), h=h(x)
%
%--------------------------------------------------------------------------
%
% L. Dede'
% EPFL, 23 May 2013
%
% Based on the Matlab NURBS toolbox by M. Spink, 
% the Octave NURBS toolbox by M. Spink, D. Claxton, C. de Falco, R. Vazquez
% and GeoPDEs: a research tool for IsoGeometric Analysis of PDEs by 
% C. de Falco, A. Reali, R. Vazquez
%--------------------------------------------------------------------------

clear all
close all
clc

%==========================================================================

% Generate Geometry
%------------------------------------
ex6_display_nurbs_surface_anular
pause( 0.5 );

% h- p- refinement
%----------------------------------------------
nurbs = nrbdegelev( nurbs, [ 0 0 ]);
new_knots1 = [ 1 : 29 ] / 30;
new_knots2 = [ 1 : 59 ] / 60;
nurbs = nrbkntins( nurbs, { new_knots1, new_knots2 } );


figure( 25 )
nrbkntplot( nurbs ); view( 0, 90 );
pause( 0.5 );

%==========================================================================

% set phyiscal DATA, boundary conditions
%----------------------------------------
mu = @( x, y )( 2e-2 + 0 * x );   % diffusion coefficient
b = @( x, y ) ( [ y; -x ] );      % advection term
rho = @( x, y )( 2e-2 + 0 * x );   % coefficient for time dependent term
f = @( x, y )( 1e2*exp(-5e2*((x-0.5).^2+(y-1.4142).^2)) );    % source term
f_time = @( t ) ( t * ( t < 1 ) + 1 * ( t >= 1 ) .* ( t <= 2 )  );  

drchlt_sides = [ 4 ];   % indexes for Dirichlet faces
nmnn_sides   = [ 1 2 3];   % indexes for Neuman faces

g = @( x, y, ind )( 0 + 0 * x );
h = @( x, y, ind )( 0 + 0 * x );    % Neuman data (for different boundaries)

drchlt_imposition_type = 'int'; % 'L2' = L2 projection, 
                                % 'int' = interpolation at control points

% time and time step
%--------------------
Tf = 4;
dt = 0.1;
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
matrix_B = op_vel_dot_gradu_v_tp( space, space, msh, b );
matrix_K = matrix_A + matrix_B;        
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
u = zeros( space.ndof, 1 );

% save vtk file
output_file_name_n = strcat( output_file_name, num2str( 10000 ) );
fprintf( 'The result is saved in the file: %s \n', output_file_name_n );
sp_to_vtk( u .* ( abs( u ) > 1e-11 ), space, geometry, vtk_pts, output_file_name_n, 'u' );


% Loop over time
%==========================================
for n = 0 : Nt - 1
    
    u_old = u;
    
    time_n = n * dt;
    time_np1 = ( n + 1 ) * dt;
    
    rhs_time_n = f_time( time_n ) * rhs_F + rhs_Neuman;
    rhs_time_np1 = f_time( time_np1 ) * rhs_F + rhs_Neuman;
    
    matrix = matrix_M + theta * dt * matrix_K;
    rhs = ( matrix_M - ( 1 - theta ) * dt * matrix_K ) * u_old ...
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
    %---------------
    output_file_name_n = strcat( output_file_name, num2str( 10000 + n + 1 ) );
    fprintf( 'The result is saved in the file: %s \n', output_file_name_n );
    sp_to_vtk( u .* ( abs( u ) > 1e-11 ), space, geometry, vtk_pts, output_file_name_n, 'u' );

    
end
              

