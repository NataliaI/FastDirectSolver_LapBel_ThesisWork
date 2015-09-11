%--------------------------------------------------------------------------
% Ex 7.2 : Plots a stream function for a multiply connected domain 
%        : Outputs Figure 7.2 
%--------------------------------------------------------------------------%--------------------------------------------------------------------------

close all; clear; clc;

% Settings for direct solver:
%----------------------------

% number of points/recursive levels :
N = 3000;
disp(['N = ' num2str(N)]);
nbox_max = 50;
h = 2*pi/N;

acc      = 1e-10; 
dist_rel = 1.6; 
nproxy   = 60; 

% Locations of islands
contour_c_ph = [pi/6; pi/8];
contour_c_th = [pi/20; 1.2];
contour_c    = [cos(contour_c_ph).*sin(contour_c_th) ...
                sin(contour_c_ph).*sin(contour_c_th) ...
                cos(contour_c_th)];

% island paramaters
flag_geom = 'star';
params = [0.2 5 0.4;
          0.2 5 0.3];
M = size(params,1);

% location of vortices and boundary data
gamma = 1;

% vortex in Omega : 
th0 = 1.2; 
ph0 = pi/3;
x0 = [cos(ph0).*sin(th0) sin(ph0).*sin(th0) cos(th0)];

% vortex on one island :  
th1 = pi/30; 
ph1 = pi/30;
x1 = [cos(ph1).*sin(th1) sin(ph1).*sin(th1) cos(th1)];

boundary = @(x,y,z) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2)/4)...
                  +(-gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2)/4);
%--------------------------------------------------------------------------


% constructing contour : 
[C,t] = OMNICONT_construct_contour_MC(flag_geom, params, contour_c, N);

% constructing second contour slightly away from boundary for plotting
[Cs,t] = OMNICONT_construct_stretched_contour_MC(flag_geom,params,contour_c,N, 3*h);

% compressing system
time_comp = tic;
[Atilde] = OMNICONT_compress_HSS_dsym_green_MC(C,nbox_max,acc,dist_rel, nproxy,N,M,contour_c);
disp(['Time for compression: ' num2str(toc(time_comp))]);

% RHS
g = boundary(C(1,:),C(2,:), C(3,:));
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
az = 150.9; el = 43.2;
v = [az el];

% Grid over sphere:
[th,ph]=meshgrid(linspace(0,pi,800),linspace(0,2*pi,1600));
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

% evaluating solution:
u = zeros(size(GridMat));
u(GridMat==0) = nan;

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

u = u + (-gamma/(2*pi))*log(sqrt((xg-x0(1)).^2+(yg-x0(2)).^2+(zg-x0(3)).^2)/4)...
       +(gamma/(2*pi))*log(sqrt((xg-x1(1)).^2+(yg-x1(2)).^2+(zg-x1(3)).^2)/4);

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
% export_fig Ex7-2TwoStarsStreamFn -m4 -png

%% CONTOURS
% Evaluating all the way up to boundary for contour plot 

[C,t] = OMNICONT_construct_contour_MC(flag_geom, params, contour_c, N);

% FIND GRIDMAT 

GridMat=zeros(size(xg));
% Evaluating integral (representation formula with sigma=1) to determine
% which grid points lie inside and outside the stretched boundary contour
for i = 1:M*N
    GridMat = GridMat+(h*C(19,i)./(2*pi)).*((xg-C(1,i)).*(C(7,i))+(yg-C(2,i)).*(C(8,i))+(zg-C(3,i)).*(C(9,i)))...
                                          ./((xg-C(1,i)).^2+(yg-C(2,i)).^2+(zg-C(3,i)).^2);
end  

GridMat(GridMat>=0)=1; % OMEGA
GridMat(GridMat<0)=0;  % ISLAND

% evaluating solution:
u = zeros(size(GridMat));
u(GridMat==0) = nan;

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

u = u + (-gamma/(2*pi))*log(sqrt((xg-x0(1)).^2+(yg-x0(2)).^2+(zg-x0(3)).^2)/4)...
       +(gamma/(2*pi))*log(sqrt((xg-x1(1)).^2+(yg-x1(2)).^2+(zg-x1(3)).^2)/4);

   
% PLOTTING CONTOURS
%------------------
figure; hold on 
xi = (xg+1i*yg)./(1-zg);

pmin = min(u(:));
%pmax = max(u(:));
pmax = 0.7;
 
pmin1 = 0.00021;
pmax1 = 0.01; 
pmin2 = 0.01;
pmax2 = 0.04;
pmin3 = 0.04;
pmax3 = 0.4; 

pinc2 = (pmax2 - pmin2)/12; 
plevels2 = pmin2:pinc2:pmax2; 
 
pinc1 = (pmax1 - pmin1)/12;
plevels1 = pmin1:pinc1:pmax1;
 
pinc3 = (pmax3 - pmin3)/30; 
plevels3 = pmin3:pinc3:pmax3;  

pinc = (pmax - pmin)/100; 
plevels = pmin:pinc:pmax; 

contour(real(xi), imag(xi), u,[plevels1 plevels2 plevels3],'LineWidth',2);

Cxi=(C(1,:)+1i*C(2,:))./(1-C(3,:));
plot(real(Cxi(1:N)), imag(Cxi(1:N)),'k','LineWidth',3);
plot(real(Cxi(N+1:end)), imag(Cxi(N+1:end)),'k','LineWidth',3);
colormap(jet); 
set(gca,'FontSize',16)
axis equal; 
xlabel('Re(\xi)');
ylabel('Im(\xi)');
set(gcf, 'Color', 'w');
axis([-10 6 -10 6]);

% export_fig Ex7-2TwoStarsStreamFnStereo -m4 -png

