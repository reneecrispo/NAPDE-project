%--------------------------------------------------------------------------
% Solve the Poisson equation in 2D by NURBS-based IGA, Galerkin method
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
display_nurbs_surface_square;
pause( 0.5 );

% h- p- refinement
%----------------------------------------------
% p=1, h-refinement
nurbs = nrbdegelev( nurbs, [ 0 0 ]);
new_knots = [ 1 : 9 ] / 10;
nurbs = nrbkntins( nurbs, { new_knots, new_knots } );

% k-refinement (starting from 1 element)
% nurbs = nrbdegelev( nurbs, [ 1 1 ]);
% new_knots = [ 1 : 9 ] / 10;
% nurbs = nrbkntins( nurbs, { new_knots, new_knots } );


figure( 25 )
nrbkntplot( nurbs ); view( 0, 90 );
pause( 0.5 );

%==========================================================================

% set phyiscal DATA, boundary conditions
%----------------------------------------
mu = @( x, y )( 1 + 0 * x );   % diffusion coefficient
f = @( x, y )( pi^2 * sin( pi * x ) .* ( 3 * ( sin( pi * y ) ).^2 - 2 * ( cos( pi * y ) ).^2 ) );    % source term

drchlt_sides = [ 1 2 3 4 ];   % indexes for Dirichlet faces
nmnn_sides   = [];   % indexes for Neuman faces

g = @( x, y, ind )( 0 + 0 * x );    % Dirichlet data (for different boundaries)
h = @( x, y, ind )( 0 + 0 * x );    % Neuman data (for different boundaries)

% set exact solution, if available
%-----------------------------------
uex = @( x, y ) ( sin( pi * x ) .* ( sin( pi * y ) ).^2 );
graduex = @( x, y ) cat( 1, ...
                   reshape( pi * cos( pi * x ) .* ( sin( pi * y ) ).^2, [1, size(x)]), ...
                   reshape( 2 * pi * sin( pi * x ) .* sin( pi * y ) .* cos( pi * y ), [1, size(x)]) );
               
% output settings
%------------------
output_file_name = 'Test2D';
n_pts_viz = 50;  % number of points for visualization

               
%==========================================================================

% Generate Space & mesh infos
%-----------------------------
geometry = geo_load( nurbs );
knots = geometry.nurbs.knots;
[ qn, qw ] = msh_set_quad_nodes( knots, ...
                                 msh_gauss_nodes( geometry.nurbs.order ) );
msh = msh_2d( knots, qn, qw, geometry );
space = sp_nurbs_2d( geometry.nurbs, msh );

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
  hval = reshape( h( x, y, iside ), msh_side.nqn, msh_side.nel );

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


% Compute errors, if exact solution is available
%------------------------------------------------
[ error_h1, error_l2 ] = sp_h1_error( space, msh, u, uex, graduex );
fprintf( 'Error in norm L2: %2.16f \n', error_l2 );
fprintf( 'Error in norm H1: %2.16f \n', error_h1 );


% prepare solution for visualization with Paraview
%---------------------------------------------------
vtk_pts = { linspace( knots{ 1 }( 1 ), knots{ 1 }( end ), n_pts_viz ), ...
            linspace( knots{ 2 }( 1 ), knots{ 2 }( end ), n_pts_viz ) };
                      
fprintf( 'The result is saved in the file: %s \n', output_file_name );
sp_to_vtk( u, space, geometry, vtk_pts, output_file_name, 'u' );




