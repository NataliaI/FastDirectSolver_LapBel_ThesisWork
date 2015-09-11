%--------------------------------------------------------------------------
% Ex 6-1 Instantaneous Point Vortex - Spherical Cap - Convergence of Error
% This script creates a convergence plot for the error in the stream
% function(Figure 6.1)
%--------------------------------------------------------------------------

close all; clear; clc;

% Settings:
%--------------------------------------------------------------------------

% Location of Cap and Point Vortex:
%----------------------------------
% phi   - azimuthal [0,2pi]
% theta - elevation [0, pi]

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
gamma = 1;       % vortex strength

flag_geom = 'cap';
params    = [sin(th_c) th_c];
contour_c = [0 0 1];

N_vec        = [2.^(3:8)];
nbox_max_vec = [2 2^2 2^3 2^3 2^4 2^4];
err          = zeros(size(N_vec));

% Direct Solver Setttings:
%--------------------------
acc      = 1e-10;
dist_rel = 1.5;
nproxy   = 50;
%--------------------------------------------------------------------------

for i = 1:length(N_vec)
    
    N = N_vec(i);
    nbox_max = nbox_max_vec(i);
    
    % constructing curve
    [C] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % boundary data and RHS
    boundary = @(x,y,z) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2))...
                      - (gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2));
    
    g = boundary(C(1,:),C(2,:),C(3,:));
    b = (2*g)';
    
    [NODES,Atilde] = OMNICONT_compress_HSS_dsym_green(C,nbox_max,acc,dist_rel,nproxy);
    
    btilde = zeros(size(Atilde,1),1);
    btilde(1:N) = b;
    
    %----umfpack--------
    [L, U, P, Q, R] = lu(Atilde);
    sigma_tilde = Q*(U\(L\(P*(R\btilde))));
    sigma = sigma_tilde(1:N);
    
    % Evaluating solution
    [th_e,ph_e] = meshgrid(linspace(5*pi/6, 5*pi/6,1),linspace(0,2*pi,100));
    xe = cos(ph_e).*sin(th_e);
    ye = sin(ph_e).*sin(th_e);
    ze = cos(th_e);
    
    h=2*pi/N;
    
    u=0;
    for j=1:N
        u=u+(h*sigma(j).*C(19,j)/(2*pi)).*(((xe-C(1,j)).*(C(7,j))+(ye-C(2,j)).*(C(8,j))+(ze-C(3,j)).*(C(9,j))))...
            ./((xe-C(1,j)).^2+(ye-C(2,j)).^2+(ze-C(3,j)).^2);
    end
    
    psi_IE= (-gamma/(2*pi))*log(sqrt((xe-x0(1)).^2+(ye-x0(2)).^2+(ze-x0(3)).^2))...
           + (gamma/(2*pi))*log(sqrt((xe-x1(1)).^2+(ye-x1(2)).^2+(ze-x1(3)).^2))...
           + u;
    
    % Exact Solution - Crowdy
    %--------------------------
    r_st=sin(th_c)/(1-cos(th_c)); % r_c in stereographic plane
    
    % stereographic projection
    z=cot(th_e/2).*exp(1i*ph_e);   % solution evaluation grid
    z0=cot(th_0/2).*exp(1i*ph_0);  % point vortex
    
    % solution on spherical cap
    zeta=1i*((r_st-z)./(r_st+z));
    zeta_0=1i*((r_st-z0)./(r_st+z0));
    
    psi_Crowdy=(-gamma/(2*pi))*log(abs((zeta-zeta_0)./(zeta-conj(zeta_0))));
    
    err(i)=norm(abs(psi_Crowdy-psi_IE), inf);
    
end

% Plot of Domain: 
%------------------
figure; hold on; 
[th_g,ph_g]=meshgrid(linspace(0, pi,50),linspace(0,2*pi,100));
xg=cos(ph_g).*sin(th_g);
yg=sin(ph_g).*sin(th_g);
zg=cos(th_g);

surf(xg, yg, zg); shading interp; alpha(0.9); axis equal; 
h1 = plot3(C(1,:),C(2,:),C(3,:)); 
h2 = plot3(xe,ye,ze);
leg = legend([h1 h2], 'boundary', 'evaluating solution','Location','SouthEast');
set(leg, 'FontSize',14); 
title('Domain');
view([86.4 6]);

% Convergence Plot: 
%-----------------------
b1=[0 0.4470 0.7410];
figure; hold on 
semilogy(N_vec, err,'Color',b1,'LineWidth',2);
semilogy(N_vec, err,'k.','MarkerSize',10);

set(gca, 'XScale', 'linear', 'YScale', 'log')
set(gca,'FontSize',14);
title('Error');
xlabel('N');
ylabel('{|| \cdot ||}_\infty');
set(gcf, 'Color', 'w');

%*** for exporting: 
% addpath('export_fig-master'); 
% export_fig Ex6-1CapError -pdf 


