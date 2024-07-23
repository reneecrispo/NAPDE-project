%--------------------------------------------------------------------------
% plot NURBS surface in R^3 with knot insertion and order elevation
%
% L. Dede'
% EPFL, 13 March 2014
%
% Based on the Matlab NURBS toolbox by M. Spink and
% the Octave NURBS toolbox by M. Spink, D. Claxton, C. de Falco, R. Vazquez
%--------------------------------------------------------------------------

clear all
close all
clc

% knots vectors in parametric direction 1 and 2
knots_1 = [ 0 0 1 1 ];
knots_2 = [ 0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1 ];
% number of basis functions in parametric direction 1 & 2
ni = 2;
nj = 9;

% specify control points and weights
xyz_p = zeros( 3, ni, nj );
ww = ones( ni, nj );
% surface properties
ttt = 0; % thick
H = 3;   % height

% First: fix parametric direction 1. Second: move along direction 2
xyz_p( 1 : 3, 1, 1 ) = [ 1; 0; 0 ]; 
ww( 1, 1 ) = 1 ;
xyz_p( 1 : 3, 1, 2 ) = [ 1; 1; 0.0 ]; 
ww( 1, 2 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 1, 3 ) = [ 0; 1; 0.0 ]; 
ww( 1, 3 ) = 1;

xyz_p( 1 : 3, 1, 4 ) = [ -1; 1; 0.0 ]; 
ww( 1, 4 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 1, 5 ) = [ -1; 0; 0.0 ]; 
ww( 1, 5 ) = 1;

xyz_p( 1 : 3, 1, 6 ) = [ -1; -1; 0.0 ]; 
ww( 1, 6 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 1, 7 ) = [ 0; -1; 0.0 ]; 
ww( 1, 7 ) = 1;

xyz_p( 1 : 3, 1, 8 ) = [ 1; -1; 0.0 ]; 
ww( 1, 8 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 1, 9 ) = [ 1; 0; 0.0 ]; 
ww( 1, 9 ) = 1;

% second row in parametric direction 1 
xyz_p( 1 : 3, 2, 1 ) = [ 1+ttt; 0; H ]; 
ww( 2, 1 ) = 1 ;
xyz_p( 1 : 3, 2, 2 ) = [ 1+ttt; 1+ttt; H ]; 
ww( 2, 2 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 2, 3 ) = [ 0; 1+ttt; H ]; 
ww( 2, 3 ) = 1;

xyz_p( 1 : 3, 2, 4 ) = [ -1-ttt; 1+ttt; H ]; 
ww( 2, 4 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 2, 5 ) = [ -1-ttt; 0; H ]; 
ww( 2, 5 ) = 1;

xyz_p( 1 : 3, 2, 6 ) = [ -1-ttt; -1-ttt; H ]; 
ww( 2, 6 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 2, 7 ) = [ 0; -1-ttt; H ]; 
ww( 2, 7 ) = 1;

xyz_p( 1 : 3, 2, 8 ) = [ 1+ttt; -1-ttt; H ]; 
ww( 2, 8 ) = sqrt( 2 ) / 2 ;
xyz_p( 1 : 3, 2, 9 ) = [ 1+ttt; 0; H ]; 
ww( 2, 9 ) = 1;


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


%==========================================================================
% knot insertion or order elevation       
%==========================================================================
% nurbs = nrbdegelev( nurbs, [ 1 1 ]);    % order elevation
% nurbs = nrbkntins( nurbs, { [ 0.4 ], [1/3 3/5] });   % knot insertion

%==========================================================================
% Do not modify below this line
%==========================================================================
%
knots_1 = nurbs.knots{1};
knots_2 = nurbs.knots{2};

ni = nurbs.number(1);
nj = nurbs.number(2);

xx = [];
yy = [];
zz = [];
for i = 1 : ni
    for j = 1 : nj
        xyz_p( 1 : 3, i, j ) = nurbs.coefs( 1 : 3, i, j ) / nurbs.coefs( 4, i, j );
            xx = [ xx, xyz_p( 1, i, j ) ];
            yy = [ yy, xyz_p( 2, i, j ) ];
            zz = [ zz, xyz_p( 3, i, j ) ];               
    end
end


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

