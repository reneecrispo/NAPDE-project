%--------------------------------------------------------------------------
% display NURBS surface in R^2, square
%
% L. Dede'
% EPFL, 11 April 2013
%
% Based on the Matlab NURBS toolbox by M. Spink and
% the Octave NURBS toolbox by M. Spink, D. Claxton, C. de Falco, R. Vazquez
%--------------------------------------------------------------------------

clear all
close all
clc

% knots vectors in parametric direction 1 and 2
knots_1 = [ 0 0 1 1 ];
knots_2 = [ 0 0 1 1 ];
% number of basis functions in parametric direction 1 & 2
ni = 2;
nj = 2;

% specify control points and weights
xyz_p = zeros( 3, ni, nj );
ww = ones( ni, nj );


% First: fix parametric direction 1. Second: move along direction 2
xyz_p( 1 : 3, 1, 1 ) = [ 0; 0; 0 ]; 
ww( 1, 1 ) = 1;
xyz_p( 1 : 3, 1, 2 ) = [ 0; 1; 0.0 ]; 
ww( 1, 2 ) = 1;

% second row in parametric direction 1 
xyz_p( 1 : 3, 2, 1 ) = [ 1; 0; 0 ]; 
ww( 2, 1 ) = 1;
xyz_p( 1 : 3, 2, 2 ) = [ 1; 1; 0 ]; 
ww( 2, 2 ) = 1;


%==========================================================================
% Do not modify below this line
%==========================================================================

% fix weights & control points
%--------------------------------
coefs = zeros( 4, ni, nj );

xx = [];
yy = [];
zz = [];
for i = 1 : ni
    for j = 1 : nj
       coefs( 1 : 3, i, j ) = xyz_p( 1 : 3, i, j ) * ww( i, j );
       coefs( 4, i, j ) = ww( i, j );
       xx = [ xx, xyz_p( 1, i, j ) ];
       yy = [ yy, xyz_p( 2, i, j ) ];
       zz = [ zz, xyz_p( 3, i, j ) ];       
    end    
end

% nurbs
%---------
knots = { knots_1, knots_2 };
nurbs = nrbmak( coefs, knots );
p = nurbs.order(1)-1;
q = nurbs.order(2)-1;

% plot
%---------
figure( 20 );
nrbplot( nurbs, [ 100 200 ], 'light', 'on' );

box on
hold on
plot3( xx, yy, zz, 'ok', 'MarkerSize', 10, 'MarkerFaceColor','r' );
for i = 1 : ni - 1
    for j = 1 : nj
        plot3( [ xyz_p( 1, i, j ) xyz_p( 1, i + 1, j ) ], ...
               [ xyz_p( 2, i, j ) xyz_p( 2, i + 1, j ) ], ...   
               [ xyz_p( 3, i, j ) xyz_p( 3, i + 1, j ) ], 'r' );
    end
end
for i = 1 : ni
    for j = 1 : nj - 1
        plot3( [ xyz_p( 1, i, j ) xyz_p( 1, i, j + 1 ) ], ...
               [ xyz_p( 2, i, j ) xyz_p( 2, i, j + 1 ) ], ...   
               [ xyz_p( 3, i, j ) xyz_p( 3, i, j + 1 ) ], 'r' );
    end
end

for i = p+1 : ni+1
  uv = linspace( knots_2( 1), knots_2( end), 401 );
  pv = nrbeval( nurbs, { knots_1( i ), uv } );
  plot3( pv(1,:), pv(2,:), pv(3,:), '-k','LineWidth',2);    
end
for j = q+1 : nj+1
  uv = linspace( knots_1( 1), knots_1( end), 401 );
  pv = nrbeval( nurbs, { uv, knots_2( j ) } );
  plot3( pv(1,:), pv(2,:), pv(3,:), '-k','LineWidth',2);    
end


hold off

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
axis equal; axis tight; 
grid


grid on
view(0,90)
axis([-0.05 1.05 -0.05 1.05])