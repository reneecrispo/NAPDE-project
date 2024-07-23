%--------------------------------------------------------------------------
% plot NURBS basis functions in the univariate parametric space
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

% basis information
%-----------------------
% Knot vector
X = [ 0 0 0 0.2 0.4 0.6 0.8 1 1 1];
% Order
p = 2;
% Number of basis functions
n = 5 + p;
% Weights
w = ones( 1, n );

%==========================================================================
% Do not modify below this line
%==========================================================================

% eval basis
%-------------------------
X = sort(X);
n_uv = 1001;
uv = linspace( X(1)+1e2*eps, X(end)-1e2*eps, n_uv );
ii = findspan(n-1,p,uv,X);
N = basisfun( ii, uv, p, X );
iN = numbasisfun( ii, uv, p, X ) + 1;

Nv = zeros( n_uv, n );
for j = 1 : n_uv
    Nv( j, iN( j, : ) ) = N( j, : );
end

Den = Nv * w';

for j = 1 : n
    Rv( :, j ) = Nv( :, j ) * w( j ) ./ Den( : );
end


% plot basis functions
%----------------------
figure( 9 );
plot( uv, Rv, 'Linewidth', 2 ); 
axis( [X(1)-0.01 X(end)+0.01 -0.01 1.01]); grid
hold on
plot(X,zeros(length(X)),'xk', 'MarkerSize',20, 'Linewidth',2);
hold off
FFSS = 20;
xSize = 8; ySize = 12;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gca,'FontSize',FFSS,'FontName','TimesNewRoman')
xlhand = get(gca,'xlabel');
set(xlhand,'FontSize',FFSS,'FontName','TimesNewRoman')
ylhand = get(gca,'ylabel');
set(ylhand,'FontSize',FFSS,'FontName','TimesNewRoman')
tlhand = get(gca,'title');
set(tlhand,'FontSize',FFSS,'FontName','TimesNewRoman')
xlabel('\xi');
ylabel('R_i(\xi)');
title({strcat('p=',num2str(p));strcat('\Xi=',mat2str(X))})
set(gca,'XTick',[ 0 0.5 1]);
set(gca,'YTick',[0 0.5 1]);
