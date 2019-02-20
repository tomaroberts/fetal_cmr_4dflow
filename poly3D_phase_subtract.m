function [Pcorr, polySurf] = poly3D_phase_subtract(P,MASK)

%% poly3D_phase_subtract   GENERATE BACKGROUND PHASE SUBTRACTED IMAGES
%
% INPUT:
%       - P             Original phase image
%       - MASK          mask covering static region in phase image
%
% OUTPUT:
%       - Pcorr         Corrected phase image
%       - polySurf      Polynomial surface
%
%

% Data Admin
[nRow,nCol,nSl,nFrame] = size( P );
M = MASK;

% Estimate Constant Background Phase
P0 = median(P(find(M)));

% Fit Polynomial to Background Phase
[xdc,ydc,zdc] = meshgrid(1:nCol,1:nRow,1:nSl);
% [xdc,ydc,zdc] = meshgrid(voxSize(2)*(1:nCol),voxSize(1)*(1:nRow),voxSize(3)*(1:nSl));
x = repmat(xdc,[1,1,nFrame]);
y = repmat(ydc,[1,1,nFrame]);
z = repmat(zdc,[1,1,nFrame]);
polyOrder = 3;
% model = polyfitn( [x(find(M)),y(find(M)),z(find(M))], angle(P(find(M))*exp(-(1i*P0))), polyOrder );
model = polyfitn( [x(find(M)),y(find(M)),z(find(M))], P(find(M)), polyOrder );

polySurf = reshape( polyvaln( model, [xdc(:),ydc(:),zdc(:)] ), [nRow,nCol,nSl] );

% Poly-corrected Phase Image
Pcorr = P-polySurf;



end % poly3D_phase_subtract(...)