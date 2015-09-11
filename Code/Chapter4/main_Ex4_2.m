%--------------------------------------------------------------------------
% Ex 4-2 : Effect of Geometry on Performance
%          Testing the fast direct solver on four different boundary
%          curves, comparing timings, skeletons, conditioning
%--------------------------------------------------------------------------
% This script corresponds to the results presented in Chapter 4, Ex 4-2.
% Creates Figure 4.4

clear; clc; close all;

% Settings for proxy method:
N_vec_P        = [2.^(10:13)];
nbox_max_vec_P = [2^6 2^6 2^7 2^7];

acc = 1e-10;
dist_rel = 1.5;
nproxy = 50;

% Location of all contours
contour_c=[0 0 1];

%--------------------------------------------------------------------------

% Storing the following info for results
% All matrices below store info in the following order:
% (:,1)  = solution error
% (:,2)  = compression time
% (:,3)  = LU decomposition time
% (:,4)  = solution time
% (:,5)  = sum of off-diagonal ranks on the coarsest level
% (:,6)  = condition number of uncompressed matrix, infinity norm

ProxyResultsGamma1 = zeros(length(N_vec_P),6);
ProxyResultsGamma2 = zeros(length(N_vec_P),6);
ProxyResultsGamma3 = zeros(length(N_vec_P),6);
ProxyResultsGamma4 = zeros(length(N_vec_P),6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour 1 -- ellipse ---

flag_geom = 'ellipse';
a0=0.8;  b0=0.4;
params = [a0 b0];

% Exact Solution
th0 = pi/20;    ph0 = pi/20;   % location of singularity
x0 = [cos(ph0)*sin(th0) sin(ph0)*sin(th0) cos(th0)];
exact = @(x,y,z) real(1./(((x+1i*y)./(1-z))-((x0(1)+1i*x0(2))./(1-x0(3)))));

% PROXY COMPRESSION
%--------------------
disp('Gamma 1');
disp('-------');
for i=1:length(N_vec_P);
    N=N_vec_P(i);
    nbox_max=nbox_max_vec_P(i);
    display(N);
    
    % constructing curve
    [C1] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C1(1,:), C1(2,:), C1(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C1,nbox_max,acc,dist_rel,nproxy);
    
    ProxyResultsGamma1(i,2) = toc(start_comp);
    ProxyResultsGamma1(i,5) = NODES_P{9,2} + NODES_P{9,3};
    
    btilde=zeros(size(Atilde_P,1),1);
    btilde(1:N)=b;
    
    %----umfpack--------
    start_lu=tic;
    [L, U, P, Q, R]=lu(Atilde_P);
    ProxyResultsGamma1(i,3) = toc(start_lu);
    
    start_solve=tic;
    sigma_tilde_P = Q*(U\(L\(P*(R\btilde))));
    ProxyResultsGamma1(i,4) = toc(start_solve);
    sigma_P=sigma_tilde_P(1:N);
    
    % Evaluating solution
    % test contour
    [theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
    x=cos(phi).*sin(theta);
    y=sin(phi).*sin(theta);
    z=cos(theta);
    
    h=2*pi/N;
    
    u_P=0;
    for j=1:N
        u_P=u_P+(h*sigma_P(j).*C1(19,j)/(2*pi)).*(((x-C1(1,j)).*(C1(7,j))+(y-C1(2,j)).*(C1(8,j))+(z-C1(3,j)).*(C1(9,j))))...
            ./((x-C1(1,j)).^2+(y-C1(2,j)).^2+(z-C1(3,j)).^2);
    end
    
    ProxyResultsGamma1(i,1) = norm(abs(u_P-exact(x,y,z)),inf);
    
    % Obtaining condition number
    disp('building matrix');
    A = OMNICONT_construct_A_diag(C1,1:N);
    
    disp('getting condition number');
    ProxyResultsGamma1(i,6) = cond(A,inf);
    
    disp('Results Gamma 1:');
    disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
    disp(ProxyResultsGamma1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour 2 -- star 1 ---

flag_geom = 'star';
a = 0.3; w = 5; c = 0.6;
params = [a w c];

% PROXY COMPRESSION
%--------------------
disp('Gamma 2');
disp('-------');
for i=1:length(N_vec_P);
    N=N_vec_P(i);
    nbox_max=nbox_max_vec_P(i);
    display(N);
    
    % constructing curve
    [C2] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C2(1,:), C2(2,:), C2(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C2,nbox_max,acc,dist_rel,nproxy);
    
    ProxyResultsGamma2(i,2) = toc(start_comp);
    ProxyResultsGamma2(i,5) = NODES_P{9,2} + NODES_P{9,3};
    
    btilde=zeros(size(Atilde_P,1),1);
    btilde(1:N)=b;
    
    %----umfpack--------
    start_lu=tic;
    [L, U, P, Q, R]=lu(Atilde_P);
    ProxyResultsGamma2(i,3) = toc(start_lu);
    
    start_solve=tic;
    sigma_tilde_P = Q*(U\(L\(P*(R\btilde))));
    ProxyResultsGamma2(i,4) = toc(start_solve);
    sigma_P=sigma_tilde_P(1:N);
    
    % Evaluating solution
    % test contour
    [theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
    x=cos(phi).*sin(theta);
    y=sin(phi).*sin(theta);
    z=cos(theta);
    
    h=2*pi/N;
    
    u_P=0;
    for j=1:N
        u_P=u_P+(h*sigma_P(j).*C2(19,j)/(2*pi)).*(((x-C2(1,j)).*(C2(7,j))+(y-C2(2,j)).*(C2(8,j))+(z-C2(3,j)).*(C2(9,j))))...
            ./((x-C2(1,j)).^2+(y-C2(2,j)).^2+(z-C2(3,j)).^2);
    end
    
    ProxyResultsGamma2(i,1) = norm(abs(u_P-exact(x,y,z)),inf);
    
    % Obtaining condition number
    disp('building matrix');
    A = OMNICONT_construct_A_diag(C2,1:N);
    
    disp('getting condition number');
    ProxyResultsGamma2(i,6) = cond(A,inf);
    
    disp('Results Gamma 2:');
    disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
    disp(ProxyResultsGamma2)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour 3 -- star 2 ---

flag_geom = 'star';
a=0.5; w=10; c=0.6;
params = [a w c];

% PROXY COMPRESSION
%--------------------
disp('Gamma 3');
disp('-------');
for i=1:length(N_vec_P);
    N=N_vec_P(i);
    nbox_max=nbox_max_vec_P(i);
    display(N);
    
    % constructing curve
    [C3] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C3(1,:), C3(2,:), C3(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C3,nbox_max,acc,dist_rel,nproxy);
    
    ProxyResultsGamma3(i,2) = toc(start_comp);
    ProxyResultsGamma3(i,5) = NODES_P{9,2} + NODES_P{9,3};
    
    btilde=zeros(size(Atilde_P,1),1);
    btilde(1:N)=b;
    
    %----umfpack--------
    start_lu=tic;
    [L, U, P, Q, R]=lu(Atilde_P);
    ProxyResultsGamma3(i,3) = toc(start_lu);
    
    start_solve=tic;
    sigma_tilde_P = Q*(U\(L\(P*(R\btilde))));
    ProxyResultsGamma3(i,4) = toc(start_solve);
    sigma_P=sigma_tilde_P(1:N);
    
    % Evaluating solution
    % test contour
    [theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
    x=cos(phi).*sin(theta);
    y=sin(phi).*sin(theta);
    z=cos(theta);
    
    h=2*pi/N;
    
    u_P=0;
    for j=1:N
        u_P=u_P+(h*sigma_P(j).*C3(19,j)/(2*pi)).*(((x-C3(1,j)).*(C3(7,j))+(y-C3(2,j)).*(C3(8,j))+(z-C3(3,j)).*(C3(9,j))))...
            ./((x-C3(1,j)).^2+(y-C3(2,j)).^2+(z-C3(3,j)).^2);
    end
    
    ProxyResultsGamma3(i,1) = norm(abs(u_P-exact(x,y,z)),inf);
    
    % Obtaining condition number
    disp('building matrix');
    A = OMNICONT_construct_A_diag(C3,1:N);
    
    disp('getting condition number');
    ProxyResultsGamma3(i,6) = cond(A,inf);
    
    disp('Results Gamma 3:');
    disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
    disp(ProxyResultsGamma3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Contour 4 - ellipse 2

flag_geom = 'ellipse';
a0=0.8; b0=0.04;
params = [a0 b0];

% PROXY COMPRESSION
%--------------------
disp('Gamma 4');
disp('-------');
for i=1:length(N_vec_P);
    N=N_vec_P(i);
    nbox_max=nbox_max_vec_P(i);
    display(N);
    
    % constructing curve
    [C4] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C4(1,:), C4(2,:), C4(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C4,nbox_max,acc,dist_rel,nproxy);
    
    ProxyResultsGamma4(i,2) = toc(start_comp);
    ProxyResultsGamma4(i,5) = NODES_P{9,2} + NODES_P{9,3};
    
    btilde=zeros(size(Atilde_P,1),1);
    btilde(1:N)=b;
    
    %----umfpack--------
    start_lu=tic;
    [L, U, P, Q, R]=lu(Atilde_P);
    ProxyResultsGamma4(i,3) = toc(start_lu);
    
    start_solve=tic;
    sigma_tilde_P = Q*(U\(L\(P*(R\btilde))));
    ProxyResultsGamma4(i,4) = toc(start_solve);
    sigma_P=sigma_tilde_P(1:N);
    
    % Evaluating solution
    % test contour
    [theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
    x=cos(phi).*sin(theta);
    y=sin(phi).*sin(theta);
    z=cos(theta);
    
    h=2*pi/N;
    
    u_P=0;
    for j=1:N
        u_P=u_P+(h*sigma_P(j).*C4(19,j)/(2*pi)).*(((x-C4(1,j)).*(C4(7,j))+(y-C4(2,j)).*(C4(8,j))+(z-C4(3,j)).*(C4(9,j))))...
            ./((x-C4(1,j)).^2+(y-C4(2,j)).^2+(z-C4(3,j)).^2);
    end
    
    ProxyResultsGamma4(i,1) = norm(abs(u_P-exact(x,y,z)),inf);
    
    % Obtaining condition number
    disp('building matrix');
    A = OMNICONT_construct_A_diag(C4,1:N);
    
    disp('getting condition number');
    ProxyResultsGamma4(i,6) = cond(A,inf);
    
    disp('Results Gamma 4:');
    disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
    disp(ProxyResultsGamma4)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Contours:

% COLORS
b1=[0 0.4470 0.7410];
r1=[0.85 0.325 0.0980];

az=-31.5; 
el=62; 

% Ellipses: Gamma 1 and Gamma 4
figure;  hold on 

[theta_a,phi_a]=meshgrid(linspace(0, pi,100),linspace(0,2*pi,200)); 
xg=cos(phi_a).*sin(theta_a);
yg=sin(phi_a).*sin(theta_a); 
zg=cos(theta_a);

surf(xg,yg,zg); shading interp;
axis equal; alpha(0.9);
axis off; 
cmap=gray(100);
colormap(cmap(50:100,:));

h7 = plot3(C1(1,:),C1(2,:),C1(3,:),'Color',r1,'LineWidth',2); 
h8 = plot3(C4(1,:),C4(2,:),C4(3,:),'Color',b1,'LineWidth',2);
set(gcf, 'Color', 'w');
leg=legend([h7 h8], '\Gamma_1','\Gamma_4');
set(leg, 'FontSize',14);
set(leg, 'Position',[0.76 0.4 0.17 0.199]);
view([az el]);

% **** for exporting figures:
% addpath(export_fig-master)
% export_fig Ex4-2Gamma1and4 -m4 -png


% Stars: Gamma 2 and Gamma 3
figure; hold on 
surf(xg,yg,zg); shading interp;
axis equal; alpha(0.9);
axis off; 
cmap=gray(100);
colormap(cmap(50:100,:)); 

h9  = plot3(C2(1,:),C2(2,:),C2(3,:),'Color',r1,'LineWidth',2);
h10 = plot3(C3(1,:),C3(2,:),C3(3,:),'Color',b1,'LineWidth',2);

set(gcf, 'Color', 'w');
leg=legend([h9 h10], '\Gamma_2','\Gamma_3');
set(leg, 'FontSize',14);
set(leg, 'Position',[0.76 0.4 0.17 0.199]);
view([az el]);

% **** for exporting figures:
% export_fig Ex4-2Gamma2and3 -m4 -png




