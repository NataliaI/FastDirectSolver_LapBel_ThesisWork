% Fast Direct Solver for the Dirichlet Laplace-Beltrami Equation 
%--------------------------------------------------------------------------
% This script implements a fast direct solver for the Laplace-Beltrami
% equation based on codes modified from : 
% http://amath.colorado.edu/faculty/martinss/2014_CBMS/codes.html
% 
% Both the brute force and proxy point methods are implemented 
% The main functions used are: 
%   OMNICONT_construct_contour
%   OMNICONT_compress_HSS_dsym_brute
%   OMNICONT_compress_HSS_dsym_green
%
%--------------------------------------------------------------------------

clear; clc; close all;

% SETTINGS: 
%--------------------------------------------------------------------------

% BOUNDARY CURVE 
%--------------------------------------------------------------------------
% 3 options: 
% - STAR
% flag_geom = 'star'; 
% params    = [a w c];   a : degree of folding/curvature
%                        w : number of arms
%                        c : size of star
%
% - ELLIPSE
% flag_geom = 'ellipse'; 
% params    = [a0 b0];   a0, b0 : major/minor axes
%
% - CAP
% flag_geom = 'cap' 
% params    = [sin(theta_c) theta_c]; theta_c : angle of elevation of cap
%                                             : theta_c in [0, pi]
%--------------------------------------------------------------------------

% Boundary 
flag_geom = 'star';  
a = 0.3; w = 5; c = 0.6; 
params = [a w c];

% Location of contour (spherical coords)
contour_c_ph=0;
contour_c_th=pi/20;
contour_c=[cos(contour_c_ph).*sin(contour_c_th) ...
           sin(contour_c_ph).*sin(contour_c_th) ...
           cos(contour_c_th)];

% number of points
N = 2^8; 

% Direct Solver Settings: 
nbox_max = 32; % number of points in one block on finest level of recursive 
               % algorithm
acc      = 1e-10; % accuracy of ID 

dist_rel = 1.5; % for proxy method - ratio which determines how proxy 
                % circle is scaled 
nproxy   = 50;  % for proxy method - number of points on proxy cap

%--------------------------------------------------------------------------

% Constructing contour: 
[C] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);
    
% Exact Solution and RHS
x0=contour_c;
exact=@(x,y,z) real(1./(((x+1i*y)./(1-z))-((x0(1)+1i*x0(2))./(1-x0(3)))));
g=exact(C(1,:), C(2,:), C(3,:));
b=(2*g)';
    
    
% BRUTE FORCE COMPRESSION
%--------------------------------------------------------------------------
    
start_comp=tic;
[NODES_BF,Atilde_BF] = OMNICONT_compress_HSS_dsym_brute(C,nbox_max,acc);
time_comp_BF=toc(start_comp);
   
btilde_BF=zeros(size(Atilde_BF,1),1);
btilde_BF(1:N)=b;

%----umfpack-----------------
start_lu=tic;
[L, U, P, Q, R]=lu(Atilde_BF);
% R : sparse diagonal 
% P : sparse
% L : lower triangular and sparse
% Q : sparse 
time_lu_BF=toc(start_lu);
  
start_solve=tic;
sigma_tilde_BF = Q*(U\(L\(P*(R\btilde_BF))));
time_solve_BF=toc(start_solve);
sigma_BF=sigma_tilde_BF(1:N);
   
    
% Evaluating solution
% (test contour is well away from boundary)
[theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
x=cos(phi).*sin(theta);
y=sin(phi).*sin(theta);
z=cos(theta);

h=2*pi/N;
    
u_BF=0;
for j=1:N
    u_BF=u_BF+(h*sigma_BF(j).*C(19,j)/(2*pi)).*(((x-C(1,j)).*(C(7,j))+(y-C(2,j)).*(C(8,j))+(z-C(3,j)).*(C(9,j))))...
                                     ./((x-C(1,j)).^2+(y-C(2,j)).^2+(z-C(3,j)).^2);
end

err_BF=norm(abs(u_BF-exact(x,y,z)),inf);

% PRINT RESULTS
disp('Brute Force Results: ');
disp('---------------------');
disp(['Total Time for Compression(s): ' num2str(time_comp_BF) 's']);
disp(['Total Time for LU Decomposition: ' num2str(time_lu_BF) 's']);
disp(['Total Time for Solution: ' num2str(time_solve_BF) 's']);
disp(' '); 
disp(['Max error over test contour: ' num2str(err_BF)]);

figure; spy(Atilde_BF); 
title('$\tilde{A}$ for Brute Force Method', 'Interpreter','latex','FontSize',14);

% PROXY COMPRESSION
%--------------------------------------------------------------------------

start_comp=tic;
[NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C,nbox_max,acc, dist_rel, nproxy);
time_comp_P=toc(start_comp);
   
btilde_P=zeros(size(Atilde_P,1),1);
btilde_P(1:N)=b;

%----umfpack-----------------
start_lu=tic;
[L, U, P, Q, R]=lu(Atilde_P);
% R : sparse diagonal 
% P : sparse
% L : lower triangular and sparse
% Q : sparse 
time_lu_P=toc(start_lu);
  
start_solve=tic;
sigma_tilde_P = Q*(U\(L\(P*(R\btilde_P))));
time_solve_P=toc(start_solve);
sigma_P=sigma_tilde_P(1:N);
   
    
% Evaluating solution
% (test contour is well away from boundary)
[theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
x=cos(phi).*sin(theta);
y=sin(phi).*sin(theta);
z=cos(theta);

h=2*pi/N;
    
u_P=0;
for j=1:N
    u_P=u_P+(h*sigma_P(j).*C(19,j)/(2*pi)).*(((x-C(1,j)).*(C(7,j))+(y-C(2,j)).*(C(8,j))+(z-C(3,j)).*(C(9,j))))...
                                     ./((x-C(1,j)).^2+(y-C(2,j)).^2+(z-C(3,j)).^2);
end

err_P=norm(abs(u_P-exact(x,y,z)),inf);

% PRINT RESULTS
disp(' ');
disp('Proxy Results: ');
disp('-------------- ');
disp(['Total Time for Compression: ' num2str(time_comp_P) 's']);
disp(['Total Time for LU Decomposition: ' num2str(time_lu_P) 's']);
disp(['Total Time for Solution: ' num2str(time_solve_P) 's']);
disp(' '); 
disp(['Max error over test contour: ' num2str(err_P)]);

figure; spy(Atilde_P); 
title('$\tilde{A}$ for Proxy Point Method', 'Interpreter','latex','FontSize',14);

%--------------------------------------------------------------------------
% Plotting boundary curve
[th_g,ph_g]=meshgrid(linspace(0,pi,50),linspace(0,2*pi,100));
xg=cos(ph_g).*sin(th_g);
yg=sin(ph_g).*sin(th_g);
zg=cos(th_g);

figure; surf(xg,yg,zg); shading interp; axis equal; alpha(0.8);
hold on 
h1 = plot3(C(1,:),C(2,:),C(3,:));
h2 = plot3(x,y,z); 
title('Boundary Curve'); 
leg = legend([h1 h2], 'boundary', 'evaluating solution','Location', 'SouthEast');
set(leg, 'FontSize',14); 

% inward normal (tangent to sphere): 
quiver3(C(1,:),C(2,:),C(3,:),C(7,:),C(8,:),C(9,:)); 
% tangent (oriented counterclockwise)
quiver3(C(1,:),C(2,:),C(3,:),C(10,:),C(11,:),C(12,:));
% e_r
quiver3(C(1,:),C(2,:),C(3,:),C(1,:),C(2,:),C(3,:));

