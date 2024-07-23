%--------------------------------------------------------------------------
% Solve the advection-diffusion equation in 2D by NURBS-based IGA, Galerkin 
% method. Based on GeoPDEs. Imposition of the Dirichlet boundary conditions 
% by means of L2 projection or "interpolation at the control points".
%
%--------------------------------------------------------------------------
%
%   - \nabla \cdot (\mu \nabla u) + b \cdot \nabla u = f,    in \Omega
%
%                               u = g,    on \Gamma_Dirichlet
%            \mu \nabla u \cdot n = h,    on \Gamma_Neumann
% 
%--------------------------------------------------------------------------
%
% L. Dede'
% EPFL, 25 April 2013
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
ex6_display_nurbs_surface_anular;
pause( 0.5 );

% h- p- refinement
%----------------------------------------------
nurbs = nrbdegelev( nurbs, [ 0 0 ]);
new_knots1 = [ 1 : 19 ] / 20;
new_knots2 = [ 1 : 49 ] / 50;
nurbs = nrbkntins( nurbs, { new_knots1, new_knots2 } );


figure( 25 )
nrbkntplot( nurbs ); view( 0, 90 );
pause( 0.5 );

%==========================================================================

% set phyiscal DATA, boundary conditions
%----------------------------------------
mu = @( x, y )( 1e-2 + 0 * x );   % diffusion coefficient
b = @( x, y ) ( [ sin( atan( y ./ ( x + 1e-12 ) ) ); - cos( atan( y ./ ( x + 1e-12 ) ) ) ] );
f = @( x, y )( 0 + 0 * x );    % source term

drchlt_sides = [ 1 2 4 ];   % indexes for Dirichlet faces
nmnn_sides   = [ 3 ];   % indexes for Neuman faces

% g = @( x, y, ind )( sqrt( x.^2 + y.^2 ) > 1.5 );
g = @( x, y, ind )( ( ind == 4 ) + 0 * x );
h = @( x, y, ind )( 0 + 0 * x );    % Neuman data (for different boundaries)

drchlt_imposition_type = 'int'; % 'L2' = L2 projection, 
                                % 'int' = interpolation at control points

% set exact solution, if available
%-----------------------------------
% uex = @( x, y ) ( 0 + 0 * x  );
% graduex = @( x, y ) cat( 1, ...
%                   reshape( 0 + 0 * x , [1, size(x)]), ...
%                   reshape( 0 + 0 * x , [1, size(x)]) );
               
% output settings
%------------------
output_file_name = 'Test2D';
n_pts_viz = 95;  % number of points for visualization

               
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
matrix = op_gradu_gradv_tp( space, space, msh, mu ) ...
         + op_vel_dot_gradu_v_tp( space, space, msh, b );

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
u( drchlt_dofs ) = u_drchlt;

int_dofs = setdiff( 1 : space.ndof, drchlt_dofs );
rhs( int_dofs ) = rhs(int_dofs) ...
                  - matrix( int_dofs, drchlt_dofs ) * u_drchlt;


% Solve the linear system (direct method)
%------------------------------------------
u( int_dofs ) = matrix( int_dofs, int_dofs ) \ rhs( int_dofs );


% Compute errors, if exact solution is available
%------------------------------------------------
% [ error_h1, error_l2 ] = sp_h1_error( space, msh, u, uex, graduex );
% fprintf( 'Error in norm L2: %2.16f \n', error_l2 );
% fprintf( 'Error in norm H1: %2.16f \n', error_h1 );


% prepare solution for visualization with Paraview
%---------------------------------------------------
vtk_pts = { linspace( knots{ 1 }( 1 ), knots{ 1 }( end ), n_pts_viz ), ...
            linspace( knots{ 2 }( 1 ), knots{ 2 }( end ), n_pts_viz ) };
                      
fprintf( 'The result is saved in the file: %s \n', output_file_name );
sp_to_vtk( u, space, geometry, vtk_pts, output_file_name, 'u' );

