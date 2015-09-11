%--------------------------------------------------------------------------
% Chapter 7 : A fast direct solver for a multiply connected domain
%           : This script solves the Laplace-Beltrami equation for a domain
%             with two star-shaped islands
%           : nbox_max is set so that the recursive levels always stay on
%             one contour at a time
%           : Produces a solution plot - Figure 7.1a)
%           : Changed colormap so figure looks different
%--------------------------------------------------------------------------

close all; clear; clc;

% Settings for direct solver:
%----------------------------

% number of points/recursive levels :
N        = 2000;
disp(['N = ' num2str(N)]);
nbox_max = 50;
h        = 2*pi/N;

acc      = 1e-10;
dist_rel = 1.6;
nproxy   = 60;

% Locations of islands
contour_c_ph = [pi/6; pi/8];
contour_c_th = [pi/6; pi/2];
contour_c    = [cos(contour_c_ph).*sin(contour_c_th) ...
                sin(contour_c_ph).*sin(contour_c_th) ...
                cos(contour_c_th)];

% island paramaters
flag_geom = 'star';
params = [0.2 5 0.4;
          0.2 5 0.3];
M = size(params,1);

% exact solution
exact = @(x,y,z) (1/2)*(real(1./(((x+1i*y)./(1-z))-(contour_c(1,1)+1i*contour_c(1,2))./(1-contour_c(1,3)))))...
                +(1/2)*(real(1./(((x+1i*y)./(1-z))-(contour_c(2,1)+1i*contour_c(2,2))./(1-contour_c(2,3)))));

%--------------------------------------------------------------------------

% constructing contour:
[C,t] = OMNICONT_construct_contour_MC(flag_geom, params, contour_c, N);

% constructing second contour slightly away from boundary for plotting
[Cs,t] = OMNICONT_construct_stretched_contour_MC(flag_geom,params,contour_c,N, 5*h);

% compressing system
time_comp = tic;
[Atilde] = OMNICONT_compress_HSS_dsym_green_MC(C,nbox_max,acc,dist_rel, nproxy,N,M,contour_c);
disp(['Time for compression: ' num2str(toc(time_comp))]);

% RHS
g = exact(C(1,:),C(2,:), C(3,:));
b = 2*g';
btilde = zeros(size(Atilde,1),1);
btilde(1:N*M) = b;

%----umfpack-----------------
time_lu = tic; 
[L, U, P, Q, R] = lu(Atilde);
disp(['Time for LU: ' num2str(toc(time_lu))]);

time_solve = tic; 
sigma_tilde = Q*(U\(L\(P*(R\btilde))));
disp(['Time for solve: ' num2str(toc(time_solve))]);

sigma = sigma_tilde(1:M*N);

% Evaluating Solution for Plotting
%-----------------------------------

% VIEW FOR PLOTS:
%-----------------
az = 118.9; el = 16.8;
v = [az el];

% Grid over sphere:
[th,ph]=meshgrid(linspace(0,pi,500),linspace(0,2*pi,1000));
xg = cos(ph).*sin(th);
yg = sin(ph).*sin(th);
zg = cos(th);

GridMat=zeros(size(xg));
% Evaluating integral (representation formula with sigma=1) to determine
% which grid points lie inside and outside the stretched boundary contour
for i=1:M*N
    GridMat=GridMat+(h*Cs(19,i)./(2*pi)).*((xg-Cs(1,i)).*(Cs(7,i))+(yg-Cs(2,i)).*(Cs(8,i))+(zg-Cs(3,i)).*(Cs(9,i)))...
        ./((xg-Cs(1,i)).^2+(yg-Cs(2,i)).^2+(zg-Cs(3,i)).^2);
end
GridMat(GridMat>=0)=1; % OMEGA
GridMat(GridMat<0)=0;  % ISLAND


figure;
surf(xg,yg,zg,GridMat);
shading flat;
title('Grid over Domain'); axis equal;
view(v);

u = zeros(size(GridMat));
u(GridMat==0) = nan;

% evaluating solution:
for k = 1:M
    A(k) = h*sigma((k-1)*N+1:k*N)'*C(19,(k-1)*N+1:k*N)';
end

for k = 1:M*N
    u = u + (h*sigma(k).*C(19,k)/(2*pi)).*(((xg-C(1,k)).*(C(7,k))+(yg-C(2,k)).*(C(8,k))+(zg-C(3,k)).*(C(9,k))))./((xg-C(1,k)).^2+(yg-C(2,k)).^2+(zg-C(3,k)).^2);
end
for k = 2:M
    u = u + A(k)*((-1/(2*pi))*(log(sqrt((xg-contour_c(1,1)).^2+(yg-contour_c(1,2)).^2+(zg-contour_c(1,3)).^2)/4))...
        +(1/(2*pi))*log(sqrt((xg-contour_c(2,1)).^2+(yg-contour_c(2,2)).^2+(zg-contour_c(2,3)).^2)/4));
end

exactu = exact(xg, yg, zg);
exactu(GridMat==0) = nan;
Error = abs(exactu(:,2:end)-u(:,2:end));
ind = ~isnan(Error);
maxError = norm(Error(ind),inf);
disp(['Maximum Error = ' num2str(maxError)]);

% PLOTTING SOLUTION
%--------------------
% changed colormap slightly from thesis
gray_colormap1   = flipud(gray(10));
gray_colormap2   = gray_colormap1(3:end,:);
grayjet_colormap = [gray_colormap2; jet(100)];

% adjusting values to colormap
cmin = min(u(:));
cmax = max(u(:));
u_new = (100)*(u - cmin)/(cmax - cmin);
u_new = u_new + 11;
u_new(GridMat==0) = 0;

figure; surf(xg,yg,zg,u_new); shading interp; colormap(grayjet_colormap);
hold on
plot3(Cs(1,:),Cs(2,:),Cs(3,:),'k','LineWidth',4);
view(v);
set(gcf, 'Color', 'w');
axis equal; axis off;

%*** for exporting:
% addpath('export_fig-master');
% export_fig Ex7-1TwoStarsSolnLapBelt -m4 -png

