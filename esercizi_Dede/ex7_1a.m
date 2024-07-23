%--------------------------------------------------------------------------
% Solve the Poisson equation in 2D by NURBS-based IGA, Galerkin method
% Based on GeoPDEs. Plot the errors in norms L2 and H1 vs. the
% characteristic size of the mesh and the total number of degree of
% freedom. Different regularities of the basis functions are considered.
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
% EPFL, 2 April 2014
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

nurbs_initial = nurbs;

% Choice of the NURBS degree, regularity of the basis and mesh size
p_vector = [ 2 : 3 ]; % must be greater or equal than 2 for this geoemtry
n_elem_vector = 2.^[ 1 : 5 ]; % in one parametric direction
h_vector = ( 3 / 4 * pi ) ./ n_elem_vector; % characteristic mesh size for the geometry

%==========================================================================

% set phyiscal DATA, boundary conditions
%----------------------------------------
mu = @( x, y )( 1 + 0 * x );   % diffusion coefficient
f = @( x, y )( -4*pi*x.*y.*cos((pi*(x.^2 + y.^2 - 1))/3) + (4*pi^2*x.^3.*y.*sin((pi*(x.^2 + y.^2 - 1))/3))/9 + (4*pi^2*x.*y.^3.*sin((pi*(x.^2 + y.^2 - 1))/3))/9);    % source term

drchlt_sides = [ 1 2 3 4 ];   % indexes for Dirichlet faces
nmnn_sides   = [  ];   % indexes for Neuman faces

g = @( x, y, ind )( 0 + 0 * x );    % Dirichlet data (for different boundaries)
h = @( x, y, ind )( 0 + 0 * x );    % Neuman data (for different boundaries)

% set exact solution, if available
%-----------------------------------
uex = @( x, y ) ( x .* y .* sin( pi / 3 * ( x.^2 + y.^2 - 1 ) ) );
graduex = @( x, y ) cat( 1, ...
                   reshape( y.*sin((pi*(x.^2 + y.^2 - 1))/3) + (2*pi*x.^2.*y.*cos((pi*(x.^2 + y.^2 - 1))/3))/3, [1, size(x)]), ...
                   reshape( x.*sin((pi*(x.^2 + y.^2 - 1))/3) + (2*pi*x.*y.^2.*cos((pi*(x.^2 + y.^2 - 1))/3))/3, [1, size(x)]) );
               
% output settings
%------------------
output_file_name = 'Test2D';
n_pts_viz = 50;  % number of points for visualization


% Loop to generate convergence figures
%---------------------------------------
for p = p_vector
    fprintf( 'Polynomial degree, p= %d \n', p );
    regularity_vector = [ 0 : p - 1 ];
    
    err_L2_matrix = [];
    err_H1_matrix = [];
    ndof_tot_matrix = [];
    
    for k = regularity_vector
        fprintf( '\t Regularity of the basis (C^k-cont.) k = %d \n', k );    
        
        ndof_tot_vector = [];
        err_L2_vect = [];
        err_H1_vect = [];
        
        for n_elem = n_elem_vector
            fprintf( '\t\t Number of elements (in one param. dir.) = %d \n', n_elem );

            % h- p- refinement
            %----------------------------------------------
            % p=1, h-refinement
            nurbs = nrbdegelev( nurbs_initial, [ ( p - 2 ) ( p - 2 ) ]);
           
            new_knots_notrepeated = [ 1 : ( n_elem - 1 ) ] / n_elem;
            new_knots = [];
            for m = p - 1 : -1 : k
                new_knots = [ new_knots, new_knots_notrepeated ];
            end
            new_knots = sort( new_knots );
            nurbs = nrbkntins( nurbs, { new_knots, new_knots } );

           
            ndof_tot_vector = [ ndof_tot_vector, nurbs.number( 1 ) * nurbs.number( 2 ) ];
            
            % figure( 250 )
            % nrbkntplot( nurbs ); view( 0, 90 );
            % pause( 0.5 );

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
            fprintf( '\t\t\t Error in norm L2: %2.16f \n', error_l2 );
            fprintf( '\t\t\t Error in norm H1: %2.16f \n', error_h1 );
            
            err_L2_vect = [ err_L2_vect, error_l2 ];
            err_H1_vect = [ err_H1_vect, error_h1 ];
            
        end    
    
        err_L2_matrix = [ err_L2_matrix; err_L2_vect ];
        err_H1_matrix = [ err_H1_matrix; err_H1_vect ];
        ndof_tot_matrix = [ ndof_tot_matrix; ndof_tot_vector ];
    
    end
    
    % plot errors, figures
    %----------------------
    
    % error L2 vs. h
    figure( 10 * p );
    
    loglog( h_vector, err_L2_matrix', '-o', ...
            h_vector, h_vector.^(p+1) * ( 3 * max( err_L2_matrix( :, end ) ) / h_vector( end )^(p+1) ), '--k', ...
            'LineWidth', 2, 'MarkerSize', 10 );
    grid on
    axis( [ 2e-2 2 1e-9 1e0 ] );
    xlabel('h');
    ylabel('err');
    M_legend = [];
    for k = regularity_vector
        M_legend = [ M_legend; strcat( 'C^', num2str( k ) ) ];
    end
    M_legend = [ M_legend; strcat( 'h^', num2str( p+1 ) ) ];
            
    legend( M_legend, 'Location','SouthEast');
    title( strcat( strcat( 'Error norm L2 vs. h, p = ', num2str( p ) ), ', basis C^k' ) );
    
    FontSize = 15;
    set( gca, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    xlhand = get( gca, 'xlabel' );
    set( xlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    ylhand = get( gca, 'ylabel' );
    set( ylhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    tlhand = get( gca, 'title' );
    set( tlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    
    % error L2 vs. ndof
    figure( 10 * p + 1 );

    loglog( ndof_tot_matrix', err_L2_matrix', '-o', ...
            'LineWidth', 2, 'MarkerSize', 10 );
    grid on
    axis( [ 1e1 1e5 1e-9 1e0 ] );
    xlabel('ndof_{tot}');
    ylabel('err');
    M_legend = [];
    for k = regularity_vector
        M_legend = [ M_legend; strcat( 'C^', num2str( k ) ) ];
    end
            
    legend( M_legend, 'Location','NorthEast');
    title( strcat( strcat( 'Error norm L2 vs. ndof, p = ', num2str( p ) ), ', basis C^k' ) );
    
    
    FontSize = 15;
    set( gca, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    xlhand = get( gca, 'xlabel' );
    set( xlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    ylhand = get( gca, 'ylabel' );
    set( ylhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    tlhand = get( gca, 'title' );
    set( tlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    
    % error H1 vs. h
    figure( 10 * p + 2 );
    
    loglog( h_vector, err_H1_matrix', '-o', ...
            h_vector, h_vector.^(p) * ( 3 * max( err_H1_matrix( :, end ) ) / h_vector( end )^(p) ), '--k', ...
            'LineWidth', 2, 'MarkerSize', 10 );
    grid on
    axis( [ 2e-2 2 1e-9 1e0 ] );
    xlabel('h');
    ylabel('err');
    M_legend = [];
    for k = regularity_vector
        M_legend = [ M_legend; strcat( 'C^', num2str( k ) ) ];
    end
    M_legend = [ M_legend; strcat( 'h^', num2str( p ) ) ];
            
    legend( M_legend, 'Location','SouthEast');
    title( strcat( strcat( 'Error norm H1 vs. h, p = ', num2str( p ) ), ', basis C^k' ) );
    
    FontSize = 15;
    set( gca, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    xlhand = get( gca, 'xlabel' );
    set( xlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    ylhand = get( gca, 'ylabel' );
    set( ylhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    tlhand = get( gca, 'title' );
    set( tlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    
    % error H1 vs. ndof
    figure( 10 * p + 3 );

    loglog( ndof_tot_matrix', err_H1_matrix', '-o', ...
            'LineWidth', 2, 'MarkerSize', 10 );
    grid on
    axis( [ 1e1 1e5 1e-9 1e0 ] );
    xlabel('ndof_{tot}');
    ylabel('err');
    M_legend = [];
    for k = regularity_vector
        M_legend = [ M_legend; strcat( 'C^', num2str( k ) ) ];
    end
            
    legend( M_legend, 'Location','NorthEast');
    title( strcat( strcat( 'Error norm H1 vs. ndof, p = ', num2str( p ) ), ', basis C^k' ) );
    
    
    FontSize = 15;
    set( gca, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    xlhand = get( gca, 'xlabel' );
    set( xlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    ylhand = get( gca, 'ylabel' );
    set( ylhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    tlhand = get( gca, 'title' );
    set( tlhand, 'FontSize', FontSize, 'FontName', 'TimesNewRoman' )
    

    pause
    
end

% prepare solution for visualization with Paraview
%---------------------------------------------------
vtk_pts = { linspace( knots{ 1 }( 1 ), knots{ 1 }( end ), n_pts_viz ), ...
            linspace( knots{ 2 }( 1 ), knots{ 2 }( end ), n_pts_viz ) };
                      
fprintf( 'The result is saved in the file: %s \n', output_file_name );
sp_to_vtk( u, space, geometry, vtk_pts, output_file_name, 'u' );




