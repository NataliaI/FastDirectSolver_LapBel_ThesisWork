%--------------------------------------------------------------------------
% Exs 6-5 and 6-7 Point Vortex Trajectories - Star-Shaped Boundary
% This script plots point vortex trajectories for a star-shaped domain
% Method 1: Constructs Contours of \hat{Psi}(x0,x0)
% - implemented first
% - outputs Figure 6.6
% Method 2: Solves IVP for trajectory
% - implemented second
% - uses spectral deferrred correction
% - outpus Figure 6.10
%--------------------------------------------------------------------------

close all; clear; clc;

% Settings:

% Location of Point Vortex on Island:
%---------------------------------------------
% phi   - azimuthal [0,2pi]
% theta - elevation [0, pi]

th_1 = pi/20; ph_1 = pi/20;   % point vortex on island
x1   = [cos(ph_1)*sin(th_1) sin(ph_1)*sin(th_1) cos(th_1)];

% Boundary: 
%-----------
gamma = 1; 

flag_geom = 'star'; 
a = 0.3; w = 5; c = 0.6; 
params    = [a w c]; 
contour_c = [0 0 1]; 

N = 2^11;
disp(['N = ' num2str(N)]);
nbox_max = 2^6;

% Direct Solver Setttings:
%--------------------------
acc      = 1e-10;
dist_rel = 1.5;
nproxy   = 50;
h        = 2*pi/N;

% Settings for Plots:
%---------------------
% colors:
gray1=[0.827451 0.827451 0.827451];
gray2=[0.662745 0.662745 0.662745];
b1=[0 0.4470 0.7410]; 

% view: 
az=-23.5; el=32;
v=[az el];
%--------------------------------------------------------------------------

% constructing curve
[C] = OMNICONT_construct_contour(flag_geom, params, contour_c, N);

% constructing second curve slightly away from boundary for plotting/evaluating
% solution 
[Cs] = OMNICONT_construct_stretched_contour(flag_geom, params, contour_c, N, 3*h);

% boundary data 
boundary = @(x,y,z,x0,x1) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2))...
                         -(gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2));

% COMPRESSING Atilde
%---------------------
time_comp=tic;
[NODES,Atilde] = OMNICONT_compress_HSS_dsym_green(C,nbox_max,acc,dist_rel,nproxy);
disp(['Time for compression: ' num2str(toc(time_comp))]);

btilde=zeros(size(Atilde,1),1);

%----umfpack---------------
[L, U, P, Q, R]=lu(Atilde);

% SOLVING FOR DIFFERENT RIGHT HAND SIDES (POINT VORTEX LOCATIONS x0)
%--------------------------------------------------------------------

% CREATING GRID OF VORTEX LOCATIONS (x0) OVER SPHERE

n = 200; % number of points on grid     
[th_0,ph_0] = meshgrid(linspace(0, pi,n),linspace(0,2*pi,2*n)); 
x0 = cos(ph_0).*sin(th_0);
y0 = sin(ph_0).*sin(th_0); 
z0 = cos(th_0);
    
GridMat=zeros(size(x0));
% Evaluating integral (representation formula with sigma=1) to determine
% which grid points lie inside and outside the bounary contour
for i=1:N
    GridMat = GridMat+(h*Cs(19,i)./(2*pi)).*((x0-Cs(1,i)).*(Cs(7,i))+(y0-Cs(2,i)).*(Cs(8,i))+(z0-Cs(3,i)).*(Cs(9,i)))...
                      ./((x0-Cs(1,i)).^2+(y0-Cs(2,i)).^2+(z0-Cs(3,i)).^2);
end
 
GridMat(GridMat>=0) = 1; % OMEGA
GridMat(GridMat<0) = 0;  % ISLAND
 
psi_hat = zeros(size(x0)); psi_hat(GridMat==0) = nan; 

time_solve = 0;
count      = 1;  % count for number of grid points 
    
for i=1:size(x0,1)
    if mod(i,10)==0
        disp(['i = ' num2str(i) ' out of ' num2str(size(x0,1))]);
    end
    for j=1:size(x0,2)
        pv_x0 = [x0(i,j) y0(i,j) z0(i,j)];
        g     = boundary(C(1,:),C(2,:), C(3,:), pv_x0, x1);
              
        b = 2*g';
        btilde(1:N) = b;
        
        start_solve=tic; 
        sigma_tilde = Q*(U\(L\(P*(R\btilde))));
        
        time_solve = time_solve + toc(start_solve);
        
        sigma = sigma_tilde(1:N);

        psi_hat(i,j) = psi_hat(i,j) + (1/2)*(sum((h*sigma'.*C(19,:)./(2*pi)).*(((pv_x0(1)-C(1,:)).*(C(7,:))+(pv_x0(2)-C(2,:)).*(C(8,:))+(pv_x0(3)-C(3,:)).*(C(9,:))))...
                      ./((pv_x0(1)-C(1,:)).^2+(pv_x0(2)-C(2,:)).^2+(pv_x0(3)-C(3,:)).^2))...
                      + (gamma/(2*pi))*log(sqrt((pv_x0(1)-x1(1)).^2+(pv_x0(2)-x1(2)).^2+(pv_x0(3)-x1(3)).^2))+(gamma/(4*pi))*log(2));
        
        % extra constant to match Crowdy's solution:
        psi_hat(i,j) = 2*pi*psi_hat(i,j)+(3/4)*log(1/2);
        count=count+1;
    end
end

disp(['Total number of grid points: ' num2str(count)]);
disp(['Total time to solve system for all grid points = ' num2str(time_solve)]);

% PLOTTING REGULAR PART OF STREAM FUNCTION AND CONTOURS
%--------------------------------------------------------
% Making two plots of each
% Second is used in spectral deferred coreection later in code

% Regular Stream Function on Sphere
gray_colormap1      = flipud(gray(10));
gray_colormap2      = gray_colormap1(3:end,:);
grayparula_colormap = [gray_colormap2; parula(100)];

% adjusting values to colormap
cmin = min(psi_hat(:));
cmax = max(psi_hat(:));
psi_hat_new = (100)*(psi_hat - cmin)/(cmax - cmin);
psi_hat_new = psi_hat_new + 11; 
psi_hat_new(GridMat==0) = 0; 

for i = 1:2
    
h1 = figure; hold on 
surf(x0,y0,z0,psi_hat_new); shading interp; colormap(grayparula_colormap);
plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',3);
axis equal; axis off; view(v); 
set(gcf, 'Color', 'w'); 

% Contours of Stream Function in Stereographic Plane
h2 = figure; hold on 
xi      = (x0+1i*y0)./(1-z0);  
pmin    = -1.9; % obtained from looking at values of psi_hat 
pmax    = max(psi_hat(:));
pinc    = ((pmax) - (pmin))/20;
plevels = pmin:pinc:pmax;

[CC, hh] = contour(real(xi),imag(xi), psi_hat, plevels, 'LineWidth', 2);
xic = (C(1,:) + 1i*C(2,:))./(1 - C(3,:));
plot(real(xic),imag(xic),'k','LineWidth',3);
axis equal; 
axis([-5 5 -5 5]);
set(gca, 'FontSize',14); 
xlabel('Re(\xi)'); 
ylabel('Im(\xi)'); 
set(gcf, 'Color', 'w');

end

%*** for exporting:
% addpath('export_fig-master');
figure(1); 
% export_fig Ex6-5StarPsiHat -m4 -png
figure(2); colorbar; 
% export_fig Ex6-5StarPsiHatStereo -pdf

%% SPECTRAL DEFERRED CORRECTION (SDC)
%--------------------------------------------------------------------------
% solving IVP for vortex trajectory with SDC (Section 5.3.2):
% dx0(t)/dt = F(t, x0(t))
%--------------------------------------------------------------------------

% Settings: 
%------------
th_0 = 13*pi/40; ph_0 = 0; % initial position of vortex
T  = 15; 
Mp = 10; 
J  = Mp-1;
npanels = 15; 
btilde = zeros(size(Atilde,1),1);

ICx0=[cos(ph_0)*sin(th_0) sin(ph_0)*sin(th_0) cos(th_0)];
    
[sol,time] = SDC(ICx0, x1, gamma, C, T, J, Mp, npanels, N, boundary, btilde, L, U, P, Q, R);
    
% PLOTTING SOLUTION in STEREOGRAPHIC PLANE for Figure 6.10a
%------------------------------------------------------------

xi   = (x0 + 1i*y0)./(1 - z0); 
soli = (sol(:,1) + 1i*sol(:,2))./(1 - sol(:,3));
    
figure(h2); colormap(gray2); 
plot(real(soli), imag(soli), 'b', 'LineWidth',2); 
plot(real(soli), imag(soli), 'b.', 'MarkerSize',9); 
% export_fig Ex6-7StarSDCStereo -pdf 

% PLOTTING SOLUTIUON ON SPHERE for Figure 6.10b
%----------------------------------------------
figure(h1); 
plot3(sol(:,1), sol(:,2), sol(:,3), 'b', 'LineWidth', 3);
% arrows were drawn manually in plot editor 
% export_fig Ex6-7StarSDC -m4 -png 



