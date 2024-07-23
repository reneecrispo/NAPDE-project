%--------------------------------------------------------------------------
% display NURBS solid
%
% L. Dede'
% EPFL, 6 March 2014
%
% Based on the Matlab NURBS toolbox by M. Spink and
% the Octave NURBS toolbox by M. Spink, D. Claxton, C. de Falco, R. Vazquez
%--------------------------------------------------------------------------

clear all
close all
clc

% knots vectors in parametric direction 1, 2, and 3
knots_1 = [ 0 0 1 1 ];
knots_2 = [ 0 0 0 1 1 1 ];
knots_3 = [ 0 0 1 1 ];

% number of basis functions in parametric directions 1, 2, and 3
ni = 2;
nj = 3;
nk = 2;

% specify control points and weights
xyz_p = zeros( 3, ni, nj, nk );
ww = ones( ni, nj, nk );

% fix \xi_1, along \xi_3 and then \xi_2
xyz_p( 1 : 3, 1, 1, 1 ) = [ 1; 0; 0 ]; 
ww( 1, 1, 1 ) = 1;
xyz_p( 1 : 3, 1, 1, 2 ) = [ 1; 0; 1 ]; 
ww( 1, 1, 2 ) = 1 ;

xyz_p( 1 : 3, 1, 2, 1 ) = [ 1; 1; 0 ]; 
ww( 1, 2, 1 ) = sqrt( 2 ) / 2;
xyz_p( 1 : 3, 1, 2, 2 ) = [ 1; 1; 1 ]; 
ww( 1, 2, 2 ) = sqrt( 2 ) / 2;

xyz_p( 1 : 3, 1, 3, 1 ) = [ 0; 1; 0 ]; 
ww( 1, 3, 1 ) = 1;
xyz_p( 1 : 3, 1, 3, 2 ) = [ 0; 1; 1 ]; 
ww( 1, 3, 2 ) = 1;


% fix \xi_1, along \xi_3 and then \xi_2
xyz_p( 1 : 3, 2, 1, 1 ) = [ 2; 0; 0 ]; 
ww( 2, 1, 1 ) = 1;
xyz_p( 1 : 3, 2, 1, 2 ) = [ 2; 0; 1 ]; 
ww( 2, 1, 2 ) = 1 ;

xyz_p( 1 : 3, 2, 2, 1 ) = [ 2; 2; 0.375 ]; 
ww( 2, 2, 1 ) = sqrt( 2 ) / 2;
xyz_p( 1 : 3, 2, 2, 2 ) = [ 2; 2; 0.625 ]; 
ww( 2, 2, 2 ) = sqrt( 2 ) / 2;

xyz_p( 1 : 3, 2, 3, 1 ) = [ 0; 2; 0 ]; 
ww( 2, 3, 1 ) = 1;
xyz_p( 1 : 3, 2, 3, 2 ) = [ 0; 2; 1 ]; 
ww( 2, 3, 2 ) = 1;


%==========================================================================
% Do not modify below this line
%==========================================================================

% fix weights & control points
%--------------------------------
coefs = zeros( 4, ni, nj, nk );

xx = [];
yy = [];
zz = [];
for i = 1 : ni
    for j = 1 : nj
        for k = 1 : nk
            coefs( 1 : 3, i, j, k ) = xyz_p( 1 : 3, i, j, k ) * ww( i, j, k );
            coefs( 4, i, j, k ) = ww( i, j, k );
            xx = [ xx, xyz_p( 1, i, j, k ) ];
            yy = [ yy, xyz_p( 2, i, j, k ) ];
            zz = [ zz, xyz_p( 3, i, j, k ) ];       
        end
    end
end    

% nurbs
%-----------
knots = { knots_1, knots_2, knots_3 };
nurbs = nrbmak( coefs, knots );

% plot
%--------
figure( 10 );
nrbplot( nurbs, [ 75 75 75 ], 'light', 'on' );

box on
hold on
plot3( xx, yy, zz, 'ok', 'MarkerSize', 10, 'MarkerFaceColor','r' );
for i = 1 : ni - 1
    for j = 1 : nj
        for k = 1 : nk
        plot3( [ xyz_p( 1, i, j, k ) xyz_p( 1, i + 1, j, k ) ], ...
               [ xyz_p( 2, i, j, k ) xyz_p( 2, i + 1, j, k ) ], ...   
               [ xyz_p( 3, i, j, k ) xyz_p( 3, i + 1, j, k ) ], 'r' );
        end
    end
end
for i = 1 : ni
    for j = 1 : nj - 1
        for k = 1 : nk
        plot3( [ xyz_p( 1, i, j, k ) xyz_p( 1, i, j + 1, k ) ], ...
               [ xyz_p( 2, i, j, k ) xyz_p( 2, i, j + 1, k ) ], ...   
               [ xyz_p( 3, i, j, k ) xyz_p( 3, i, j + 1, k ) ], 'r' );
        end
    end
end
for i = 1 : ni
    for j = 1 : nj
        for k = 1 : nk - 1
        plot3( [ xyz_p( 1, i, j, k ) xyz_p( 1, i, j, k + 1 ) ], ...
               [ xyz_p( 2, i, j, k ) xyz_p( 2, i, j, k + 1 ) ], ...   
               [ xyz_p( 3, i, j, k ) xyz_p( 3, i, j, k + 1 ) ], 'r' );
        end
    end
end


p = nurbs.order(1)-1;
q = nurbs.order(2)-1;
r = nurbs.order(3)-1;
%
uv = linspace( knots_2( 1), knots_2( end), 401 );
for i = p+1 : ni+1  
  for k = r+1 : nk+1        
        pv = nrbeval( nurbs, { knots_1( i ), uv, knots_3( k ) } );
        plot3( pv(1,:), pv(2,:), pv(3,:), '-k','LineWidth',2);    
  end
end
%
uv = linspace( knots_1( 1), knots_1( end), 401 );
for j = q+1 : nj+1
  for k = r+1 : nk+1        
        pv = nrbeval( nurbs, { uv, knots_2( j ), knots_3( k ) } );
        plot3( pv(1,:), pv(2,:), pv(3,:), '-k','LineWidth',2);    
  end
end
%
uv = linspace( knots_3( 1), knots_3( end), 401 );
for i = p+1 : ni+1
  for j = q+1 : nj+1        
        pv = nrbeval( nurbs, { knots_1( i ), knots_2( j ), uv } );
        plot3( pv(1,:), pv(2,:), pv(3,:), '-k','LineWidth',2);    
  end
end


hold off

axis equal
colormap gray

FFSS = 20;
xSize = 8; ySize = 12;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gca,'FontSize',FFSS,'FontName','TimesNewRoman')
xlhand = get(gca,'xlabel');
set(xlhand,'FontSize',FFSS,'FontName','TimesNewRoman')
ylhand = get(gca,'ylabel');
set(ylhand,'FontSize',FFSS,'FontName','TimesNewRoman')
zlhand = get(gca,'zlabel');
set(zlhand,'FontSize',FFSS,'FontName','TimesNewRoman')

tlhand = get(gca,'title');
set(tlhand,'FontSize',FFSS,'FontName','TimesNewRoman')
xlabel('x');
ylabel('y');
zlabel('z');

grid on
axis([-0.2 2.2 -0.2 2.2 -0.2 1.2])