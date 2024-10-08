%--------------------------------------------------------------------------
% plot NURBS curve in R^3
%
% L. Dede'
% EPFL, 26 March 2014
%
% Based on the Matlab NURBS toolbox by M. Spink and
% the Octave NURBS toolbox by M. Spink, D. Claxton, C. de Falco, R. Vazquez
%--------------------------------------------------------------------------

clear all
close all
clc

% knot vector
knots = [ 0, 0, 0, 0.125, 0.125, 0.25, 0.25, 0.375, 0.375, 0.5, 0.5, 0.625, 0.625, 0.75, 0.75, 0.875, 0.875, 1, 1, 1 ];
% number of basis functions
ni = 17;

% control points and weights
xyz_p = zeros( 3, ni, 1 );
ww = ones( ni, 1 );

xyz_p( 1 : 3, 1, 1 ) = [ 1; 0; 0 ]; 
ww( 1, 1) = 1;

xyz_p( 1 : 3, 2, 1 ) = [ 1; 1; 0.25 ]; 
ww( 2, 1 ) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 3, 1 ) = [ 0; 1; 0.5 ]; 
ww( 3, 1 ) = 1;

xyz_p( 1 : 3, 4, 1 ) = [ -1; 1; 0.75 ]; 
ww( 4, 1) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 5, 1 ) = [ -1; 0; 1 ]; 
ww( 5, 1 ) = 1;

xyz_p( 1 : 3, 6, 1 ) = [ -1; -1; 1.25 ]; 
ww( 6, 1 ) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 7, 1 ) = [ 0; -1; 1.5 ]; 
ww( 7, 1 ) = 1;

xyz_p( 1 : 3, 8, 1 ) = [ 1; -1; 1.75 ]; 
ww( 8, 1 ) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 9, 1 ) = [ 1; 0; 2 ]; 
ww( 9, 1) = 1;

xyz_p( 1 : 3, 10, 1 ) = [ 1; 1; 2.25 ]; 
ww( 10, 1 ) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 11, 1 ) = [ 0; 1; 2.5 ]; 
ww( 11, 1 ) = 1;

xyz_p( 1 : 3, 12, 1 ) = [ -1; 1; 2.75 ]; 
ww( 12, 1) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 13, 1 ) = [ -1; 0; 3 ]; 
ww( 13, 1 ) = 1;

xyz_p( 1 : 3, 14, 1 ) = [ -1; -1; 3.25 ]; 
ww( 14, 1 ) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 15, 1 ) = [ 0; -1; 3.5 ]; 
ww( 15, 1 ) = 1;

xyz_p( 1 : 3, 16, 1 ) = [ 1; -1; 3.75 ]; 
ww( 16, 1 ) = 1 / sqrt( 2 );

xyz_p( 1 : 3, 17, 1 ) = [ 1; 0; 4 ]; 
ww( 17, 1) = 1;

%==========================================================================
% Do not modify below this line
%==========================================================================

% fix weights & control points
%--------------------------------
coefs = zeros( 4, ni, 1 );

xx = [];
yy = [];
zz = [];
for i = 1 : ni
       coefs( 1 : 3, i, 1 ) = xyz_p( 1 : 3, i, 1 ) * ww( i, 1 );
       coefs( 4, i, 1 ) = ww( i, 1 );
       xx = [ xx, xyz_p( 1, i, 1 ) ];
       yy = [ yy, xyz_p( 2, i, 1 ) ];
       zz = [ zz, xyz_p( 3, i, 1 ) ];       
end    

% nurbs
%-----------
nurbs = nrbmak( coefs, knots );

% plot
%--------
figure( 10 );
nrbplot( nurbs, 1001, 'light', 'on' );
set(findobj(gca,'Type','line'),'LineWidth',2)

box on
hold on
plot3( xx, yy, zz, 'or', 'MarkerSize', 10, 'MarkerFaceColor','r' );
for i = 1 : ni - 1
    plot3( [ xyz_p( 1, i, 1 ) xyz_p( 1, i + 1, 1 ) ], ...
           [ xyz_p( 2, i, 1 ) xyz_p( 2, i + 1, 1 ) ], ...   
           [ xyz_p( 3, i, 1 ) xyz_p( 3, i + 1, 1 ) ], 'r' );
end
pknts = nrbeval( nurbs, knots );
plot3( pknts(1,:), pknts(2,:), pknts(3,:), 'xk', 'MarkerSize', 20, ...
       'LineWidth', 2, 'MarkerFaceColor','k' );

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
grid on
