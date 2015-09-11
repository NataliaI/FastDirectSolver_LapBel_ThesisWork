%--------------------------------------------------------------------------
% Ex 6-2 Instantaneous Point Vortex - Star-Shaped Domain - Solution Plots
% This script plots the instantaneous stream function for point vortex
% motion in a star-shaped domain
% Ouputs Figure 6.3
%--------------------------------------------------------------------------

close all; clear; clc; 

% Settings: 
%--------------------------------------------------------------------------

% Location of Cap and Point Vortex:
%----------------------------------
% phi   - azimuthal [0,2pi]
% theta - elevation [0, pi]

th_0=pi/3;  ph_0=3*pi/2;   % point vortex in Omega
th_1=pi/8;  ph_1=3*pi/2;   % point vortex on island

% point vortex in Omega (Cartesian coords):
x0=[cos(ph_0)*sin(th_0) sin(ph_0)*sin(th_0) cos(th_0)];

% point vortex on cap (Cartesian coords):
x1=[cos(ph_1)*sin(th_1) sin(ph_1)*sin(th_1) cos(th_1)];

% Boundary:
%----------
gamma = 1; % vortex strength

flag_geom = 'star'; 
a = 0.3; w = 5; c = 0.6; 
params = [a w c]; 
contour_c = [0 0 1]; 

N=2^11; 
disp(['N = ' num2str(N)]);
nbox_max = 2^7; 

% Direct Solver Setttings:
%--------------------------
acc      = 1e-10;
dist_rel = 1.5;
nproxy   = 50;
h        = 2*pi/N; 

%--------------------------------------------------------------------------

% constructing curve
[C]  = OMNICONT_construct_contour(flag_geom, params, contour_c, N);
% constructing second curve slightly away from boundary for plotting
[Cs] = OMNICONT_construct_stretched_contour(flag_geom, params, contour_c, N, 5*h);

% boundary data and RHS
boundary = @(x,y,z) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2)/4)...
             - (gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2)/4);
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

% Evaluating Solution for Plotting
%-----------------------------------

% VIEW FOR PLOTS: 
%-----------------
az=-23.5;  
el=32;
v=[az el];

% Grid over sphere: 
[th,ph]=meshgrid(linspace(0,pi,500),linspace(0,2*pi,1000)); 
xg=cos(ph).*sin(th); 
yg=sin(ph).*sin(th);
zg=cos(th);

h = 2*pi/N;

GridMat=zeros(size(xg));
% Evaluating integral (representation formula with sigma=1) to determine
% which grid points lie inside and outside the stretched boundary contour
for i=1:N
    GridMat=GridMat+(h*Cs(19,i)./(2*pi)).*((xg-Cs(1,i)).*(Cs(7,i))+(yg-Cs(2,i)).*(Cs(8,i))+(zg-Cs(3,i)).*(Cs(9,i)))...
                     ./((xg-Cs(1,i)).^2+(yg-Cs(2,i)).^2+(zg-Cs(3,i)).^2);
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
psi_IE_new=(100)*(psi_IE-cmin)/(cmax-cmin);
psi_IE_new=psi_IE_new+11; 
psi_IE_new(GridMat==0)=0; 

figure; surf(xg,yg,zg,psi_IE_new); shading interp; 
colormap(grayjet_colormap); 
hold on 
plot3(Cs(1,:),Cs(2,:),Cs(3,:),'k','LineWidth',4);
view(v); 
axis equal; 
axis off;
set(gcf, 'Color', 'w'); 
%*** for exporting: 
% addpath('export_fig-master'); 
% export_fig Ex6-2StarStreamFn -m4 -png 
 
% PLOTTING Stream Function Psi(x,x0) in Stereographic Plane 
%------------------------------------------------------------

% Evaluating all the way up to the boundary: 
GridMat=zeros(size(xg));
% Evaluating integral (representation formula with sigma=1) to determine
% which grid points lie inside and outside the bounary contour
for i=1:N
    GridMat=GridMat+(h*C(19,i)./(2*pi)).*((xg-C(1,i)).*(C(7,i))+(yg-C(2,i)).*(C(8,i))+(zg-C(3,i)).*(C(9,i)))...
                     ./((xg-C(1,i)).^2+(yg-C(2,i)).^2+(zg-C(3,i)).^2);
end
GridMat(GridMat>=0)=1; % OMEGA
GridMat(GridMat<0)=0;  % ISLAND

u=zeros(size(GridMat));
u(GridMat==0)=nan; 

for j=1:N
    u=u+(h*sigma(j).*C(19,j)/(2*pi)).*(((xg-C(1,j)).*(C(7,j))+(yg-C(2,j)).*(C(8,j))+(zg-C(3,j)).*(C(9,j))))...
        ./((xg-C(1,j)).^2+(yg-C(2,j)).^2+(zg-C(3,j)).^2);
end

psi_IE = (-gamma/(2*pi))*log(sqrt((xg-x0(1)).^2+(yg-x0(2)).^2+(zg-x0(3)).^2)/4)...
         +(gamma/(2*pi))*log(sqrt((xg-x1(1)).^2+(yg-x1(2)).^2+(zg-x1(3)).^2)/4)...
         + u;

xi=(xg+1i*yg)./(1-zg); 
figure; hold on
pmin=1e-4; pmax=0.7; pmed=0.08; pmed1=0.01; pmed2=0.05; % determined from values 
                                                        % of stream function 

pinc1=(pmed1-pmin)/6;
pinc2=(pmed2-pmed1)/6;
pinc3=(pmed-pmed2)/3; 
pinc4=(pmax-pmed)/30;
plevels1=pmin:pinc1:pmed1; 
plevels2=pmed1:pinc2:pmed2; 
plevels3=pmed2:pinc3:pmed; 
plevels4=pmed:pinc4:pmax; 

contour(real(xi), imag(xi), psi_IE,[plevels1 plevels2 plevels3 plevels4],'LineWidth',2);
Cxi=(C(1,:)+1i*C(2,:))./(1-C(3,:));

plot(real(Cxi), imag(Cxi),'k','LineWidth',3);
colormap(jet); 
axis equal; axis([-5 5 -5 5]);
set(gcf, 'Color', 'w'); 
set(gca,'FontSize',16)
xlabel('Re(\xi)');
ylabel('Im(\xi)');

% export_fig Ex6-2StarStreamFnStereo -pdf 

