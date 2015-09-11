%--------------------------------------------------------------------------
% Ex 6-1 Instantaneous Point Vortex - Spherical Cap - Solution Plots
% This script plots the instantaneous stream function for point vortex
% motion in a spherical cap domain
% Ouputs Figure 6.2
%--------------------------------------------------------------------------

close all; clear; clc; 

% Settings: 
%--------------------------------------------------------------------------

% Location of Cap and Point Vortex:
%----------------------------------
% phi   - azimuthal [0,2pi]
% theta - elevation [0, pi]

% Defining boundary and location of vortices
%------------------------------------------------

th_c = pi/6;                   % spherical cap
th_0 = pi/3;  ph_0 = 3*pi/2;   % point vortex in Omega
th_1 = pi/8;  ph_1 = 3*pi/2;   % point vortex on cap

% point vortex in Omega (Cartesian coords):
x0 = [cos(ph_0)*sin(th_0) sin(ph_0)*sin(th_0) cos(th_0)];

% point vortex on cap (Cartesian coords):
x1 = [cos(ph_1)*sin(th_1) sin(ph_1)*sin(th_1) cos(th_1)];

% Boundary:
%----------
r_c = sin(th_c);  % radius of cap
gamma = 1;        % vortex strength

flag_geom = 'cap';
params    = [sin(th_c) th_c];
contour_c = [0 0 1];

N        = 700;
disp(['N = ' num2str(N)]);
nbox_max = 25;

% Direct Solver Setttings:
%--------------------------
acc      = 1e-10;
dist_rel = 1.5;
nproxy   = 50;
%--------------------------------------------------------------------------

% constructing curve
[C] = OMNICONT_construct_contour(flag_geom, params, contour_c, N);

% boundary data and RHS
boundary = @(x,y,z) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2))...
                  - (gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2));
    
g = boundary(C(1,:),C(2,:),C(3,:));
b = (2*g)';

time_comp=tic; 
[NODES,Atilde] = OMNICONT_compress_HSS_dsym_green(C,nbox_max,acc,dist_rel,nproxy);
disp(['Time for compression: ' num2str(toc(time_comp))]); 

btilde=zeros(size(Atilde,1),1);
btilde(1:N)=b;

%----umfpack--------
[L, U, P, Q, R]=lu(Atilde);
sigma_tilde = Q*(U\(L\(P*(R\btilde))));
sigma=sigma_tilde(1:N);
    
% Evaluating solution (just to test error)
[th_e,ph_e]=meshgrid(linspace(5*pi/6, 5*pi/6,1),linspace(0,2*pi,100)); 
xe=cos(ph_e).*sin(th_e); 
ye=sin(ph_e).*sin(th_e);
ze=cos(th_e);

h=2*pi/N; 

u=0;
for j=1:N
    u=u+(h*sigma(j).*C(19,j)/(2*pi)).*(((xe-C(1,j)).*(C(7,j))+(ye-C(2,j)).*(C(8,j))+(ze-C(3,j)).*(C(9,j))))./((xe-C(1,j)).^2+(ye-C(2,j)).^2+(ze-C(3,j)).^2);
end
  
psi_IE=(-gamma/(2*pi))*log(sqrt((xe-x0(1)).^2+(ye-x0(2)).^2+(ze-x0(3)).^2)/4)...
       +(gamma/(2*pi))*log(sqrt((xe-x1(1)).^2+(ye-x1(2)).^2+(ze-x1(3)).^2)/4)...
       + u;

% Exact Solution - Crowdy
%--------------------------
r_st=sin(th_c)/(1-cos(th_c)); % r_c in stereographic plane

% stereographic projection
xi=cot(th_e/2).*exp(1i*ph_e);   % solution evaluation grid
xi_0=cot(th_0/2).*exp(1i*ph_0);  % point vortex

% solution on spherical cap
zeta=1i*((r_st-xi)./(r_st+xi));
zeta_0=1i*((r_st-xi_0)./(r_st+xi_0));

psi_Crowdy=(-gamma/(2*pi))*log(abs((zeta-zeta_0)./(zeta-conj(zeta_0))));

err=norm(abs(psi_Crowdy-psi_IE), inf);
disp(['Max error on test contour: ' num2str(err)]);

% Evaluating Solution for Plotting
%-----------------------------------

% VIEW FOR PLOTS: 
%-----------------
az=-23.5;  
el=32;
v=[az el];

% Grid over sphere: 
[th,ph]=meshgrid(linspace(0,pi,100),linspace(0,2*pi,200)); 
xg=cos(ph).*sin(th); 
yg=sin(ph).*sin(th);
zg=cos(th);

GridMat=zeros(size(xg));
% Evaluating integral (representation formula with sigma=1) to determine
% which grid points lie inside and outside the bounary contour
% For spherical cap - evaluating solution right up to boundary 
for i=1:N
    GridMat=GridMat+(h*C(19,i)./(2*pi)).*((xg-C(1,i)).*(C(7,i))+(yg-C(2,i)).*(C(8,i))+(zg-C(3,i)).*(C(9,i)))...
                     ./((xg-C(1,i)).^2+(yg-C(2,i)).^2+(zg-C(3,i)).^2);
end
GridMat(GridMat>=0)=1; % OMEGA
GridMat(GridMat<0)=0;  % ISLAND

figure; 
surf(xg,yg,zg,GridMat); 
shading flat; 
title('Grid over Domain'); axis equal; 
  
u=zeros(size(GridMat));
u(GridMat==0)=nan; 

for j=1:N
    u=u+(h*sigma(j).*C(19,j)/(2*pi)).*(((xg-C(1,j)).*(C(7,j))+(yg-C(2,j)).*(C(8,j))+(zg-C(3,j)).*(C(9,j))))...
        ./((xg-C(1,j)).^2+(yg-C(2,j)).^2+(zg-C(3,j)).^2);
end
  
psi_IE = (-gamma/(2*pi))*log(sqrt((xg-x0(1)).^2+(yg-x0(2)).^2+(zg-x0(3)).^2)/4)...
         +(gamma/(2*pi))*log(sqrt((xg-x1(1)).^2+(yg-x1(2)).^2+(zg-x1(3)).^2)/4)...
         + u;

% PLOTTING Stream Function Psi(x,x0)
%-------------------------------------
gray_colormap1   = flipud(gray(10));
gray_colormap2   = gray_colormap1(3:end,:);
grayjet_colormap = [gray_colormap2; jet(100)];

% adjusting values to colormap 
cmin = min(psi_IE(:));
cmax = max(psi_IE(:));
psi_IE_new = (100)*(psi_IE-cmin)/(cmax-cmin);
psi_IE_new = psi_IE_new+11; 
psi_IE_new(GridMat==0) = 0; 

figure; surf(xg,yg,zg,psi_IE_new); shading interp; 
colormap(grayjet_colormap); 
hold on 
plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',4);
view(v); 
axis equal; 
axis off;
set(gcf, 'Color', 'w'); 
%*** for exporting: 
% addpath('export_fig-master'); 
% export_fig Ex6-1CapStreamFn -m4 -png 

% PLOTTING Stream Function Psi(x,x0) in Stereographic Plane
%-----------------------------------------------------------

xi=(xg+1i*yg)./(1-zg); 
figure; hold on
pmin=0; pmax=0.8; pmed=0.2; % determined from values of stream function
pinc=(pmed-pmin)/10; 
plevels=pmin:pinc:pmed; 
pinc2=(pmax-pmed)/20;
plevels2=pmed:pinc2:pmax; 

contour(real(xi), imag(xi), psi_IE,[plevels plevels2],'LineWidth',2);
Cxi=(C(1,:)+1i*C(2,:))./(1-C(3,:));

plot(real(Cxi), imag(Cxi),'k','LineWidth',3);
colormap(jet); 
axis equal; axis([-4 4 -4 4]);
set(gcf, 'Color', 'w');
set(gca,'FontSize',16);
xlabel('Re(\xi)');
ylabel('Im(\xi)');
% export_fig Ex6-1CapStreamFnStereo -pdf 


% PLOTTING G(x,x0) (Singular Part)
%--------------------------------------
Gxx0=(-gamma/(2*pi))*log(sqrt((xg-x0(1)).^2+(yg-x0(2)).^2+(zg-x0(3)).^2)/4);

% adjusting values to colormap 
cmin = min(Gxx0(:));
cmax = max(Gxx0(:));
Gxx0_new=(100)*(Gxx0-cmin)/(cmax-cmin);
Gxx0_new=Gxx0_new+11; 
Gxx0_new(GridMat==0)=0; 

figure; surf(xg,yg,zg,Gxx0_new); shading interp; 
colormap(grayjet_colormap); 
hold on 
plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',4);
view(v); 
axis equal; 
axis off;
set(gcf, 'Color', 'w'); 
% export_fig Ex6-1CapSingPart -m4 -png 


% PLOTTING Psi_hat(x,x0) (Regular Part)
%---------------------------------------
psi_hat = (gamma/(2*pi))*log(sqrt((xg-x1(1)).^2+(yg-x1(2)).^2+(zg-x1(3)).^2)/4) + u;

% adjusting values to colormap 
cmin = min(psi_hat(:));
cmax = max(psi_hat(:));
psi_hat_new=(100)*(psi_hat-cmin)/(cmax-cmin);
psi_hat_new=psi_hat_new+11; 
psi_hat_new(GridMat==0)=0; 

figure; surf(xg,yg,zg,psi_hat_new); shading interp; 
colormap(grayjet_colormap); 
hold on 
plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',4);
view(v); 
axis equal; 
axis off;
set(gcf, 'Color', 'w'); 
% export_fig Ex6-1CapRegPart -m4 -png 
