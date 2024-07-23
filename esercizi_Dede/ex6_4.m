%--------------------------------------------------------------------------
% Solve the Poisson equation in 3D by NURBS-based IGA, Galerkin method
% Based on GeoPDEs
%
%--------------------------------------------------------------------------
%
%   - \nabla \cdot (\mu \nabla u) = f,    in \Omega
%
%                               u = g,    on \Gamma_Dirichlet
%            \mu \nabla u \cdot n = h,    on \Gamma_Neumann
% 
%--------------------------------------------------------------------------
%
% L. Dede'
% EPFL, 27 March 2014
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
display_nurbs_solid_cube;
pause( 0.5 );

% h- p- refinement
%----------------------------------------------
nurbs = nrbdegelev( nurbs, [ 1 1 1 ]);
new_knots = [ 1 : 9 ] / 10;
nurbs = nrbkntins( nurbs, { new_knots, new_knots, new_knots } );

figure( 25 )
nrbkntplot( nurbs );
pause( 0.5 );

%==========================================================================

% set phyiscal DATA, boundary conditions
%----------------------------------------
mu = @( x, y, z )( 1 + 0 * x );   % diffusion coefficient
f = @( x, y, z )( 3 * pi^2 * sin( pi * x ) .* sin( pi * y ) .* sin( pi * z ) );    % source term

drchlt_sides = [ 1 2 3 4 5 6];   % indexes for Dirichlet faces
nmnn_sides   = [];   % indexes for Neuman faces

g = @( x, y, z, ind )( 0 + 0 * x );    % Dirichlet data (for the whole boundary)
h = @( x, y, z, ind )( 0 + 0 * x );% Neuman data (for different boundaries)

% set exact solution, if available
%-----------------------------------
uex = @( x, y, z ) ( sin( pi * x ) .* sin( pi * y ) .* sin( pi * z ) );
graduex = @( x, y, z ) cat( 1, ...
                   reshape( pi * cos( pi * x ) .* sin( pi * y ) .* sin( pi * z ), [1, size(x)]), ...
                   reshape( pi * sin( pi * x ) .* cos( pi * y ) .* sin( pi * z ), [1, size(x)]), ...
                   reshape( pi * sin( pi * x ) .* sin( pi * y ) .* cos( pi * z ), [1, size(x)]) );
               
% output settings
%------------------
output_file_name = 'Test3D';
n_pts_viz = 50;  % number of points for visualization

               
%==========================================================================

% Generate Space & mesh infos
%-----------------------------
geometry = geo_load( nurbs );
knots = geometry.nurbs.knots;
[ qn, qw ] = msh_set_quad_nodes( knots, ...
                                 msh_gauss_nodes( geometry.nurbs.order ) );
msh = msh_3d( knots, qn, qw, geometry );
space = sp_nurbs_3d( geometry.nurbs, msh );

% assemble the matrix
%---------------------
matrix = op_gradu_gradv_tp( space, space, msh, mu );

% assemble rhs
%---------------
rhs = op_f_v_tp( space, msh, f );

% apply Neuman BCs
%------------------
for iside = nmnn_sides
  msh_side = msh_eval_boundary_side( msh, iside );
  sp_side  = sp_eval_boundary_side( space, msh_side );

  x = squeeze( msh_side.geo_map( 1, :, : ) );
  y = squeeze( msh_side.geo_map( 2, :, : ) );
  z = squeeze( msh_side.geo_map( 3, :, : ) );
  hval = reshape( h( x, y, z, iside ), msh_side.nqn, msh_side.nel );

  rhs( sp_side.dofs ) = rhs( sp_side.dofs ) ...
                        + op_f_v( sp_side, msh_side, hval );
end

% Apply Dirichlet BCs
%----------------------
u = zeros( space.ndof, 1 );
[ u_drchlt, drchlt_dofs ] = sp_drchlt_l2_proj( space, msh, ...
                                               g, drchlt_sides );
u( drchlt_dofs ) = u_drchlt;

int_dofs = setdiff( 1 : space.ndof, drchlt_dofs );
rhs( int_dofs ) = rhs(int_dofs) ...
                  - matrix( int_dofs, drchlt_dofs ) * u_drchlt;

% Solve the linear system (direct method)
%------------------------------------------
u( int_dofs ) = matrix( int_dofs, int_dofs ) \ rhs( int_dofs );


%==========================================================================

% Compute errors, if exact solution is available
%------------------------------------------------
[ error_h1, error_l2 ] = sp_h1_error( space, msh, u, uex, graduex );
fprintf( 'Error in norm L2: %2.10f \n', error_l2 );
fprintf( 'Error in norm H1: %2.10f \n', error_h1 );


% prepare solution for visualization with Paraview
%---------------------------------------------------
vtk_pts = { linspace( knots{ 1 }( 1 ), knots{ 1 }( end ), n_pts_viz ), ...
            linspace( knots{ 2 }( 1 ), knots{ 2 }( end ), n_pts_viz ), ...
            linspace( knots{ 3 }( 1 ), knots{ 3 }( end ), n_pts_viz ) };
                      
fprintf( 'The result is saved in the file %s \n \n', output_file_name );
sp_to_vtk( u, space, geometry, vtk_pts, output_file_name, 'u' );




