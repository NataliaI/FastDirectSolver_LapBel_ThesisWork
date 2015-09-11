
function [] = main_Ex6_4_FDSvsFMM()
%--------------------------------------------------------------------------
% Ex 6-4 : Comparing performance of Direct solver and FMM for a spherical
%          cap domain
%        : FMM routine uses FMM code in LapBel_matlab and fmmlib2d
%            - Calls the functions:
%              FMMSolve.m
%              island_geometry.m
%              build_system.m
%        : evaluating contours of regular part of stream function to
%          compare performance
%        : Outputs Figure 6.5
%--------------------------------------------------------------------------

close all; clear; clc;
addpath('matlab');  % for fmmlib mex files

az=-23.5;
el=32;
v=[az el];

b1=[0 0.4470 0.7410];
r1=[0.85 0.325 0.0980];
% view:

% Settings:
%--------------------------------------------------------------------------

% Direct Solver:
%----------------

% Location of Cap and Point Vortex:
%----------------------------------
% phi   - azimuthal [0,2pi]
% theta - elevation [0, pi]

% Defining boundary and location of vortices
%--------------------------------------------

th_c = pi/6;                     % spherical cap
th_1 = pi/20;    ph_1 = pi/20;   % point vortex on cap

% location of point vortex on cap
x1 = [cos(ph_1)*sin(th_1) sin(ph_1)*sin(th_1) cos(th_1)];

% Boundary:
%----------
gamma = 1;        % vortex strength

flag_geom = 'cap';
params    = [sin(th_c) th_c];
contour_c_ph = 0;
contour_c_th = pi/20;
contour_c    = [cos(contour_c_ph).*sin(contour_c_th) ...
    sin(contour_c_ph).*sin(contour_c_th) ...
    cos(contour_c_th)];

N        = 2^10;
nbox_max = 2^6;

acc      = 1e-10;
dist_rel = 1.5;
nproxy   = 50;
h   = 2*pi/N;

% constructing curve for fast direct solver
[C] = OMNICONT_construct_contour(flag_geom, params, contour_c, N);
% constructing curve slighly away from boundary for plotting
[Cs] = OMNICONT_construct_stretched_contour(flag_geom,params,contour_c, N, 5*h);

% boundary data for fast direct solver
boundary = @(x,y,z,x0,x1) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2))...
                         -(gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2));
% BUILD SYSTEM FOR DIRECT SOLVER
%--------------------------------
start_comp = tic;
[NODES, Atilde] = OMNICONT_compress_HSS_dsym_green(C,nbox_max,acc,dist_rel,nproxy);
time_comp =toc(start_comp); 

btilde = zeros(size(Atilde,1),1);

%----umfpack--------
start_lu = tic; 
[L, U, P, Q, R] = lu(Atilde);
time_lu = toc(start_lu); 

% SETTINGS FOR FMM:
%-------------------
% Using same variable names/notation as in FMM code
th_k = 0;  phi_k = pi/2-pi/20;

% BUILD SYSTEM for FMM:
%-----------------------
% Using parts of code from lapBel_driver.m and the function
% Laplace_Beltrami.m
% Building system for input into GMRES later
%--------------------------------------------------------------------------

% Taken from lapBel_driver.m :
%---------------------------
A  = sin(th_c);
B  = -sin(th_c);
Np = N; % number of points per island 
        % here we have one so Np = N

nbod = 1;
island_geo      = zeros(4, nbod);
island_geo(1,:) = A(1:nbod);
island_geo(2,:) = B(1:nbod);
island_geo(3,:) = th_k(1:nbod);
island_geo(4,:) = phi_k(1:nbod);

[~, nbod] = size(island_geo);
nbk = nbod * Np;

% from Laplace_Beltrami.m
%-------------------------
% Define closed curve, parametrized by alpha
[dth, alph, R_FMM, T, Normal, dsda, diagK, Ck] ...
    = island_geometry(island_geo, nbod, Np);

% stereographic projection of points
zeta_k = cart_to_zeta(Ck(1,:), Ck(2,:), Ck(3,:));
zeta   = cart_to_zeta(R_FMM(1,:), R_FMM(2,:), R_FMM(3,:));
dzeta  = (T(1,:) + 1i*T(2,:) + zeta.*T(3,:)).*dsda./(1 - R_FMM(3,:));

%
% Construct System:  [K | E ] sigma
%                    [F | D ] A_k   = f
start_build_systemFMM = tic;
[K, E, F, D] = build_system(nbod, Np, nbk, dth, R_FMM, Normal, dsda, diagK, Ck);
time_build_systemFMM = toc(start_build_systemFMM); 

%--------------------------------------------------------------------------
% Plotting contours of regular part of stream function over a fixed grid
% on the sphere
% Comparing the error between the FMM and Direct Solver
%--------------------------------------------------------------------------

n_fix = 25; % number of points on grid

[th_0,ph_0] = meshgrid(linspace(0, pi,n_fix),linspace(0,2*pi,2*n_fix));
x0 = cos(ph_0).*sin(th_0);
y0 = sin(ph_0).*sin(th_0);
z0 = cos(th_0);

GridMat = zeros(size(x0)); 
% Evaluating integral (representation formula with sigma=1) to determine
% which grid points lie inside and outside the bounary contour
for i=1:N
    GridMat = GridMat+(h*Cs(19,i)./(2*pi)).*((x0-Cs(1,i)).*(Cs(7,i))+(y0-Cs(2,i)).*(Cs(8,i))+(z0-Cs(3,i)).*(Cs(9,i)))...
                      ./((x0-Cs(1,i)).^2+(y0-Cs(2,i)).^2+(z0-Cs(3,i)).^2);
end

GridMat(GridMat>=0) = 1; % OMEGA
GridMat(GridMat<0) = 0;  % ISLAND

psi_hat_FDS = zeros(size(x0)); % psi_hat_FDS(GridMat==0) = nan; 
psi_hat_FMM = zeros(size(x0)); % psi_hat_FMM(GridMat==0) = nan; 
diff        = zeros(size(x0)); % diff(GridMat==0) =  nan; 

time_solve_FDS = 0;
time_solve_FMM = 0;

for i=1:size(x0,1)
    if mod(i,10)==0
        disp(['i = ' num2str(i) ' out of ' num2str(size(x0,1))]);
    end
    for j=1:size(x0,2)
        if GridMat(i,j)==0
            psi_hat_FDS(i,j)=nan; 
            psi_hat_FMM(i,j)=nan; 
            diff(i,j)=nan; 
        else
        pv_x0 = [x0(i,j) y0(i,j) z0(i,j)];
        g = boundary(C(1,:),C(2,:), C(3,:), pv_x0,x1);
        
        b = 2*g';
        btilde(1:N) = b;
        
        start_solve_FDS = tic;
        sigma_tilde     = Q*(U\(L\(P*(R\btilde))));
        end_solve_FDS   = toc(start_solve_FDS);
        time_solve_FDS  = time_solve_FDS + end_solve_FDS;
        sigma_FDS = sigma_tilde(1:N);
        
        % EVALUATING SOLUTION AT x0
        % GETTING SOLUTION FROM FAST DIRECT SOLVER
        psi_hat_FDS(i,j) = (1/2)*(sum((h*sigma_FDS'.*C(19,:)./(2*pi)).*(((pv_x0(1)-C(1,:)).*(C(7,:))+(pv_x0(2)-C(2,:)).*(C(8,:))+(pv_x0(3)-C(3,:)).*(C(9,:))))...
            ./((pv_x0(1)-C(1,:)).^2+(pv_x0(2)-C(2,:)).^2+(pv_x0(3)-C(3,:)).^2))...
            +(gamma/(2*pi))*log(sqrt((pv_x0(1)-x1(1)).^2+(pv_x0(2)-x1(2)).^2+(pv_x0(3)-x1(3)).^2))...
            +(gamma/(4*pi))*log(2));
        
        % GETTING SOLUTION FROM FMM
        % FMMSolve is a function adapted from Laplace_Beltrami.m and its
        % subfunctions
        [t_fmm, sigma_FMM, psi_hat_FMM(i,j)] = FMMSolve(R_FMM, zeta, nbk, D,...
            F, E, diagK, dth, dzeta,...
            Ck, Normal, dsda, acc,...
            pv_x0, x1, gamma);
        end_solve_FMM = t_fmm; 
        time_solve_FMM = time_solve_FMM + t_fmm;
        
        diff(i,j) = abs(psi_hat_FDS(i,j)-psi_hat_FMM(i,j));
        
        end
    end
end

max_diff = max(diff(:));

% PLOTTING SOLUTIONS TO psi hat FROM FDS AND FMM:

% FDS:
gray_colormap1 = flipud(gray(10));
gray_colormap2 = gray_colormap1(3:end,:);
grayparula_colormap = [gray_colormap2; parula(100)];

% adjusting values to colormap
cmin = min(psi_hat_FDS(:));
cmax = max(psi_hat_FDS(:));
psi_hat_FDS_new = (100)*(psi_hat_FDS - cmin)/(cmax - cmin);
psi_hat_FDS_new = psi_hat_FDS_new + 11; 
psi_hat_FDS_new(GridMat==0) = 0; 

figure; hold on;
surf(x0,y0,z0,psi_hat_FDS_new); shading interp; colormap(grayparula_colormap);
plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth', 4);
view(v); axis equal; axis off;
set(gcf, 'Color','w');
title('psi hat FDS');

% FMM:
% adjusting values to colormap
cmin = min(psi_hat_FMM(:));
cmax = max(psi_hat_FMM(:));
psi_hat_FMM_new = (100)*(psi_hat_FMM - cmin)/(cmax - cmin);
psi_hat_FMM_new = psi_hat_FMM_new + 11; 
psi_hat_FMM_new(GridMat==0) = 0; 

figure; hold on;
surf(x0,y0,z0,psi_hat_FMM_new); shading interp; colormap(grayparula_colormap);
plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth', 4);
view(v); axis equal; axis off;
set(gcf, 'Color','w');
title('psi hat FMM');


% TIME PLOT vs Grid Size - Figure 6.5
%-------------------------------------
% Timings for FMM and FDS for increasing grid sizes 

%nvec=[10:5:40];
nvec = [10:5:30];

time_solve_vec_FDS = zeros(size(nvec));
time_solve_vec_FMM = zeros(size(nvec));

ncount = 1;

for n = nvec
    disp(['Grid Size : ' num2str(n) ' by ' num2str(2*n)]);
    
    [th_0,ph_0] = meshgrid(linspace(th_c + 2*h, pi,n),linspace(0,2*pi,2*n));
    x0 = cos(ph_0).*sin(th_0);
    y0 = sin(ph_0).*sin(th_0);
    z0 = cos(th_0);
    
    psi_hat_FDS = zeros(size(x0));
    psi_hat_FMM = zeros(size(x0));
    
    time_solve_FDS = 0;
    time_solve_FMM = 0;
    
    for i=1:size(x0,1)
        if mod(i,10)==0
        disp(['i = ' num2str(i) ' out of ' num2str(size(x0,1))]);
        end
        for j=1:size(x0,2)
            pv_x0 = [x0(i,j) y0(i,j) z0(i,j)];
            g = boundary(C(1,:),C(2,:), C(3,:), pv_x0,x1);
        
            b = 2*g';
            btilde(1:N) = b;
            
            start_solve_FDS = tic;
            sigma_tilde     = Q*(U\(L\(P*(R\btilde))));
            time_solve_FDS  = time_solve_FDS + toc(start_solve_FDS);
            sigma_FDS       = sigma_tilde(1:N);
            
            % EVALUATING SOLUTION AT x0
            % GETTING SOLUTION FROM FAST DIRECT SOLVER
            psi_hat_FDS(i,j) = (1/2)*(sum((h*sigma_FDS'.*C(19,:)./(2*pi)).*(((pv_x0(1)-C(1,:)).*(C(7,:))+(pv_x0(2)-C(2,:)).*(C(8,:))+(pv_x0(3)-C(3,:)).*(C(9,:))))...
                                ./((pv_x0(1)-C(1,:)).^2+(pv_x0(2)-C(2,:)).^2+(pv_x0(3)-C(3,:)).^2))...
                                +(gamma/(2*pi))*log(sqrt((pv_x0(1)-x1(1)).^2+(pv_x0(2)-x1(2)).^2+(pv_x0(3)-x1(3)).^2))...
                                +(gamma/(4*pi))*log(2));
                
            % GETTING SOLUTION FROM FMM
            % FMMSolve is a function adapted from Laplace_Beltrami.m and its
            % subfunctions
            [t_fmm, sigma_FMM, psi_hat_FMM(i,j)] = FMMSolve(R_FMM, zeta, nbk, D,...
                                                            F, E, diagK, dth, dzeta,...
                                                            Ck, Normal, dsda, acc,...
                                                            pv_x0, x1, gamma);
        
            time_solve_FMM = time_solve_FMM + t_fmm;
        end
    end

    time_solve_vec_FDS(ncount) = time_solve_FDS;
    time_solve_vec_FMM(ncount) = time_solve_FMM;
    ncount = ncount+1;
end


figure; hold on
nvec_tot = nvec.*(2*nvec);
h1  = loglog(nvec_tot, time_solve_vec_FDS,'Color',b1,'LineWidth',2);
      loglog(nvec_tot, time_solve_vec_FDS,'k.','MarkerSize', 9);
h2  = loglog(nvec_tot, time_solve_vec_FMM,'Color',r1,'LineWidth',2);
      loglog(nvec_tot, time_solve_vec_FMM,'k.','MarkerSize', 9);
leg = legend([h1 h2], 'Direct Solver', 'FMM','Location','SouthEast');

set(leg, 'FontSize',14); 
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'FontSize',14);
xlabel('N_g');
ylabel('Time (s)');
set(gcf, 'Color','w'); 

%*** for exporting:
% addpath('export_fig-master');
% export_fig Ex6-4FMMVsDirect -pdf


% SUMMARY 
%---------
disp('SUMMARY'); 
disp(' '); 

disp(['N = ' num2str(N)]);

disp(['Time for compression: ' num2str(time_comp)]);
disp(['Time for LU factorization: ' num2str(time_lu)]); 
disp(['Time for building system for FMM: ' num2str(time_build_systemFMM)]);

disp('Comparing FMM and Direct Solver over fixed grid: ');
disp(['Grid Size : ' num2str(n_fix) ' by ' num2str(2*n_fix)]);
disp('Time for a Single Solve: ');
disp(['FDS : Time to apply umfpack factors : ' num2str(end_solve_FDS)]);
disp(['FMM : Time to apply GMRES to system: ' num2str(end_solve_FMM)]);
disp(['Maximum difference between solutions to FMM and FDS over fixed grid: ' num2str(max_diff)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zeta = cart_to_zeta(x,y,z)
zeta = (x + 1i*y)./(1 - z);
end

