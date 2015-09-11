%--------------------------------------------------------------------------
% Exs 6-3 and 6-6 Point Vortex Trajectories - Spherical Cap - Solution Plots
% This script plots point vortex trajectories for a spherical cap domain
% Method 1: Constructs Contours of \hat{Psi}(x0,x0)
% - implemented first
% - outputs Figure 6.4
% Method 2: Solves IVP for trajectory
% - implemented second
% - uses spectral deferrred correction
% - outpus Figure 6.7
%--------------------------------------------------------------------------

close all; clear; clc;

% Settings :
%--------------------------------------------------------------------------

% Location of Cap and Point Vortex on Cap:
%----------------------------------------
% phi   - azimuthal [0,2pi]
% theta - elevation [0, pi]

th_c = pi/6;                % spherical cap
th_1 = pi/20; ph_1=pi/20;   % point vortex on cap

% point vortex on cap (Cartesian coords):
x1 = [cos(ph_1)*sin(th_1) sin(ph_1)*sin(th_1) cos(th_1)];

% Boundary:
%----------
r_c = sin(th_c);  % radius of cap
gamma = 1;        % vortex strength

flag_geom = 'cap';
params    = [sin(th_c) th_c];
contour_c = [0 0 1];

N = 2^11;
disp(['N = ' num2str(N)]);
nbox_max = 2^6;

% Direct Solver Setttings:
%--------------------------
acc      = 1e-10;
dist_rel = 1.5;
nproxy   = 50;

% Settings for Plots:
%---------------------
% colors:
gray1=[0.827451 0.827451 0.827451];
gray2=[0.662745 0.662745 0.662745];
b1=[0 0.4470 0.7410];
r1=[0.85 0.325 0.0980];

% view:
az=-23.5; el=32;
v=[az el];
%--------------------------------------------------------------------------

% constructing curve
[C] = OMNICONT_construct_contour(flag_geom, params, contour_c, N);

% boundary data 
boundary = @(x,y,z,x0,x1) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2))...
    -(gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2));

% COMPRESSING Atilde
%---------------------
time_comp=tic;
[NODES,Atilde] = OMNICONT_compress_HSS_dsym_green(C,nbox_max,acc,dist_rel,nproxy);
disp(['Time for compression: ' num2str(toc(time_comp))]);

btilde=zeros(size(Atilde,1),1);

%----umfpack--------
[L, U, P, Q, R]=lu(Atilde);

% SOLVING FOR DIFFERENT RIGHT HAND SIDES (POINT VORTEX LOCATIONS x0)
%--------------------------------------------------------------------

% CREATING GRID OF VORTEX LOCATIONS (x0) OVER SPHERE
n = 50; % number of points on grid
h = 2*pi/N;
[th_0,ph_0] = meshgrid(linspace(th_c + 2*h, pi,n),linspace(0,2*pi,2*n));
x0 = cos(ph_0).*sin(th_0);
y0 = sin(ph_0).*sin(th_0);
z0 = cos(th_0);

psi_hat      = zeros(size(x0));
psi_hat_Crow = zeros(size(x0));
err          = zeros(size(x0));

time_solve   = 0;

disp(['Grid Size : ' num2str(n) ' by ' num2str(2*n)]);

for i=1:size(x0,1)
    if mod(i,10)==0
        disp(['i = ' num2str(i) ' out of ' num2str(size(x0,1))]);
    end
    for j=1:size(x0,2)
        pv_x0 = [x0(i,j) y0(i,j) z0(i,j)];
        g = boundary(C(1,:),C(2,:), C(3,:), pv_x0,x1);
        
        b = 2*g';
        btilde(1:N) = b;
        
        start_solve=tic;
        sigma_tilde = Q*(U\(L\(P*(R\btilde))));
        
        time_solve = time_solve + toc(start_solve);
        
        sigma = sigma_tilde(1:N);
        
        psi_hat(i,j) = psi_hat(i,j) + (1/2)*(sum((h*sigma'.*C(19,:)./(2*pi)).*(((pv_x0(1)-C(1,:)).*(C(7,:))+(pv_x0(2)-C(2,:)).*(C(8,:))+(pv_x0(3)-C(3,:)).*(C(9,:))))...
            ./((pv_x0(1)-C(1,:)).^2+(pv_x0(2)-C(2,:)).^2+(pv_x0(3)-C(3,:)).^2))...
            +(gamma/(2*pi))*log(sqrt((pv_x0(1)-x1(1)).^2+(pv_x0(2)-x1(2)).^2+(pv_x0(3)-x1(3)).^2))...
            +(gamma/(4*pi))*log(2));
        
        % extra constant to match Crowdy's solution:
        psi_hat(i,j)=2*pi*psi_hat(i,j)+(3/4)*log(1/2);
        
        % Exact Solution - Crowdy
        r_st = sin(th_c)/(1-cos(th_c));
        xi_0 = (pv_x0(1)+1i*pv_x0(2))/(1-pv_x0(3));
        zeta = 1i*(r_st-xi_0)./(r_st+xi_0);
        zetap = -2*1i*r_st./(xi_0+r_st).^2;
        
        psi_hat_Crow(i,j)=psi_hat_Crow(i,j)-(1/(2))*log(abs((1+xi_0.*conj(xi_0)).*(zetap)./(zeta-conj(zeta))));
        
        err(i,j)=err(i,j)+abs(psi_hat_Crow(i,j)-psi_hat(i,j));
    end
end

disp(['Total time to solve system for all grid points x0 = ' num2str(time_solve)]);

% PLOTTING REGULAR PART OF STREAM FUNCTION AND CONTOURS
%--------------------------------------------------------
% Making two plots of each
% Second is used in spectral deferred coreection later in code

% Regular Stream Function on Sphere
for i = 1:2
    h1 = figure; hold on
    % grid for cap:
    [th_cap,ph_cap]=meshgrid(linspace(0, th_c+2*h,25),linspace(0,2*pi,50));
    x_cap=cos(ph_cap).*sin(th_cap);
    y_cap=sin(ph_cap).*sin(th_cap);
    z_cap=cos(th_cap);
    
    surf(x0,y0,z0,psi_hat); shading interp; colormap(parula);
    freezeColors;
    surf(x_cap,y_cap,z_cap);
    colormap(gray1); shading interp;
    
    plot3(C(1,:),C(2,:),C(3,:), 'k', 'LineWidth',4);
    view(v); axis equal;  axis off;
    set(gcf, 'Color', 'w');
    
    
    % Contours of Stream Function in Stereographic Plane
    h2 = figure; hold on
    xi = (x0+1i*y0)./(1-z0);
    pmin = min(psi_hat(:));
    pmax = max(psi_hat(:));
    pinc = ((pmax)-(pmin))/20;
    plevels = pmin:pinc:pmax;
    
    [CC, hh] = contour(real(xi),imag(xi),psi_hat,plevels,'LineWidth',2);
    xic=(C(1,:)+1i*C(2,:))./(1-C(3,:));
    plot(real(xic),imag(xic),'k','LineWidth',3);
    set(gca, 'FontSize',14);
    xlabel('Re(\xi)');
    ylabel('Im(\xi)');
    set(gcf, 'Color', 'w');
    axis equal
    % export_fig Ex6-3CapPsiHatStereo -m4 -pdf
end

%*** for exporting:
% addpath('export_fig-master');
figure(1); 
% export_fig Ex6-3CapPsiHat -m4 -png
figure(2); colorbar;  
% export_fig Ex6-3CapPsiHatStereo -pdf

%% SPECTRAL DEFERRED CORRECTION (SDC)
%--------------------------------------------------------------------------
% solving IVP for vortex trajectory with SDC (Section 5.3.2):
% dx0(t)/dt = F(t, x0(t))
%--------------------------------------------------------------------------

% Settings: 
%------------
th_0 = 13*pi/40; ph_0 = 0; % initial possition of vortex
T  = 30;  % final time t=T
Mp = 10;  % number of points per panel
J  = Mp-1; % number of correction steps 
npanels = 3;
btilde = zeros(size(Atilde,1),1);

ICx0 = [cos(ph_0)*sin(th_0) sin(ph_0)*sin(th_0) cos(th_0)];

% Solving IVP: 
%-------------
[sol, ~, sol_eul] = SDC(ICx0, x1, gamma, C, T, J, Mp, npanels, N, boundary, btilde, L, U, P, Q, R); 

% PLOTTING SOLUTION in STEREOGRAPHIC PLANE for Figure 6.7a
%------------------------------------------------------------

% Stereographic Projection of Forward Euler and SDC Solutions 
xi = (x0+1i*y0)./(1-z0);

soli     = (sol(:,1) + 1i*sol(:,2))./(1 - sol(:,3));
sol_euli = (sol_eul(:,1) + 1i*sol_eul(:,2))./(1 - sol_eul(:,3));

figure(h2); colormap(gray2);
plot(real(soli), imag(soli), 'b', 'LineWidth', 2);
plot(real(soli), imag(soli), 'b.', 'MarkerSize', 9);
plot(real(sol_euli), imag(sol_euli), 'r', 'LineWidth', 2);
plot(real(sol_euli), imag(sol_euli), 'r.','MarkerSize', 9);
% export_fig Ex6-6CapSDCStereo -pdf 

% PLOTTING SOLUTIUON ON SPHERE for Figure 6.7b
%----------------------------------------------
% In order to get a smoother looking solution curve on the sphere I take
% more panels for the plot
npanels = 5;

% Solving IVP: 
%-------------
[sol, time, sol_eul] = SDC(ICx0, x1, gamma, C, T, J, Mp, npanels, N, boundary, btilde, L, U, P, Q, R); 

figure(h1); hold on
plot3(sol(:,1), sol(:,2), sol(:,3), 'Color', 'b', 'LineWidth',3);
view([-147.5078, 35.1025]);
% arrows were drawn manually in plot editor 
% export_fig Ex6-6CapSDC -m4 -png 





