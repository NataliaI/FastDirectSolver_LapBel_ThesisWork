%--------------------------------------------------------------------------
% Ex 4-1 : Brute Force Vs Proxy Point Method
%          Testing the fast direct solver for three different boundary
%          curves
%--------------------------------------------------------------------------
% This script corresponds to the results presented in Chapter 4, Ex 4-1.
% Creates Figures 4.1-4.3

clear; clc; close all;

% Settings for brute and proxy methods:

% number of points/recursive levels for brute force method:
N_vec_BF=[2.^(6:16)];
nbox_max_vec_BF=[2^4 2^4 2^5 2^5 2^6 2^6 2^7 2^7 2^7 2^7 2^7];

% number of points/recursive levels for proxy point method:
N_vec_P=[2.^(6:17)];
nbox_max_vec_P=[2^4 2^4 2^5 2^5 2^6 2^6 2^7 2^7 2^7 2^7 2^7 2^7];

acc      = 1e-10; % accuracy of ID

dist_rel = 1.5; % for proxy method - ratio which determines how proxy
                % circle is scaled
nproxy   = 50;  % for proxy method - number of points on proxy cap

% Location of all contours
contour_c=[0 0 1];

%--------------------------------------------------------------------------

% Storing the following info for results/plotting
% All matrices below store info in the following order:
% (:,1)        = solution error
% (:,2)        = compression time
% (:,3)        = LU decomposition time
% (:,4)        = solution time
% (:,5), (:,6) = off-diagonal ranks on the coarsest level

BruteResultsGamma1 = zeros(length(N_vec_BF),6);
ProxyResultsGamma1 = zeros(length(N_vec_P),6);

BruteResultsGamma2 = zeros(length(N_vec_BF),6);
ProxyResultsGamma2 = zeros(length(N_vec_P),6);

BruteResultsGamma3 = zeros(length(N_vec_BF),6);
ProxyResultsGamma3 = zeros(length(N_vec_P),6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamma 1 -- ellipse ---

flag_geom = 'ellipse';
a0 = 0.8; b0 = 0.4;
params = [a0 b0];

% Exact Solution
th0 = pi/20;    ph0 = pi/20;   % location of singularity
x0 = [cos(ph0)*sin(th0) sin(ph0)*sin(th0) cos(th0)];
exact = @(x,y,z) real(1./(((x+1i*y)./(1-z))-((x0(1)+1i*x0(2))./(1-x0(3)))));

% BRUTE FORCE COMPRESSION
%----------------------------

disp('Gamma 1: ');
disp('----------');
disp('Brute Force Compression: ');
disp(' ');

for i=1:length(N_vec_BF);
    
    N=N_vec_BF(i);
    nbox_max=nbox_max_vec_BF(i);
    disp(N);
    
    % constructing curve
    [C1] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C1(1,:), C1(2,:), C1(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_BF,Atilde_BF] = OMNICONT_compress_HSS_dsym_brute(C1,nbox_max,acc);
    
    BruteResultsGamma1(i,2) = toc(start_comp);
    BruteResultsGamma1(i,5) = NODES_BF{9,2};
    BruteResultsGamma1(i,6) = NODES_BF{9,3};
    
    btilde=zeros(size(Atilde_BF,1),1);
    btilde(1:N)=b;
    
    %----umfpack--------
    start_lu = tic;
    [L, U, P, Q, R] = lu(Atilde_BF);
    BruteResultsGamma1(i,3) = toc(start_lu);
    
    start_solve=tic;
    sigma_tilde_BF = Q*(U\(L\(P*(R\btilde))));
    BruteResultsGamma1(i,4) = toc(start_solve);
    sigma_BF=sigma_tilde_BF(1:N);
    
    % Evaluating solution
    % test contour
    [theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
    x=cos(phi).*sin(theta);
    y=sin(phi).*sin(theta);
    z=cos(theta);
    
    h=2*pi/N;
    
    u_BF=0;
    for j=1:N
        u_BF=u_BF+(h*sigma_BF(j).*C1(19,j)/(2*pi)).*(((x-C1(1,j)).*(C1(7,j))+(y-C1(2,j)).*(C1(8,j))+(z-C1(3,j)).*(C1(9,j))))...
            ./((x-C1(1,j)).^2+(y-C1(2,j)).^2+(z-C1(3,j)).^2);
    end
    
    BruteResultsGamma1(i,1) = norm(abs(u_BF-exact(x,y,z)),inf);
end

% PROXY COMPRESSION
%----------------------------

disp('Proxy Compression: ');
disp(' ');

for i=1:length(N_vec_P);
    N=N_vec_P(i);
    nbox_max=nbox_max_vec_P(i);
    disp(N);
    
    % constructing curve
    [C1] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C1(1,:), C1(2,:), C1(3,:));
    b=(2*g)';
    
    start_comp = tic;
    [NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C1,nbox_max,acc,dist_rel,nproxy);
    
    ProxyResultsGamma1(i,2) = toc(start_comp);
    ProxyResultsGamma1(i,5) = NODES_P{9,2};
    ProxyResultsGamma1(i,6) = NODES_P{9,3};
    
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour 2 -- star1 ---

flag_geom='star';
a=0.3; w=5; c=0.6;
params=[a w c];

% BRUTE FORCE COMPRESSION
%----------------------------

disp('Gamma 2: ');
disp('----------');
disp('Brute Force Compression: ');
disp(' ');

for i=1:length(N_vec_BF);
    
    N=N_vec_BF(i);
    nbox_max=nbox_max_vec_BF(i);
    disp(N);
    
    % constructing curve
    [C2] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    %RHS
    g=exact(C2(1,:), C2(2,:), C2(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_BF,Atilde_BF] = OMNICONT_compress_HSS_dsym_brute(C2,nbox_max,acc);
    
    BruteResultsGamma2(i,2) = toc(start_comp);
    BruteResultsGamma2(i,5) = NODES_BF{9,2};
    BruteResultsGamma2(i,6) = NODES_BF{9,3};
    
    btilde=zeros(size(Atilde_BF,1),1);
    btilde(1:N)=b;
    
    %----umfpack--------
    start_lu=tic;
    [L, U, P, Q, R]=lu(Atilde_BF);
    BruteResultsGamma2(i,3) = toc(start_lu);
    
    start_solve=tic;
    sigma_tilde_BF = Q*(U\(L\(P*(R\btilde))));
    BruteResultsGamma2(i,4) = toc(start_solve);
    sigma_BF=sigma_tilde_BF(1:N);
    
    % Evaluating solution
    % test contour
    [theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
    x=cos(phi).*sin(theta);
    y=sin(phi).*sin(theta);
    z=cos(theta);
    
    h=2*pi/N;
    
    u_BF=0;
    for j=1:N
        u_BF=u_BF+(h*sigma_BF(j).*C2(19,j)/(2*pi)).*(((x-C2(1,j)).*(C2(7,j))+(y-C2(2,j)).*(C2(8,j))+(z-C2(3,j)).*(C2(9,j))))...
            ./((x-C2(1,j)).^2+(y-C2(2,j)).^2+(z-C2(3,j)).^2);
    end
    
    BruteResultsGamma2(i,1)=norm(abs(u_BF-exact(x,y,z)),inf);
end

% PROXY COMPRESSION
%----------------------------

disp('Proxy Compression: ');
disp(' ');

for i=1:length(N_vec_P);
    N=N_vec_P(i);
    nbox_max=nbox_max_vec_P(i);
    disp(N);
    
    % constructing curve:
    [C2] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);
    
    % RHS
    g=exact(C2(1,:), C2(2,:), C2(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C2,nbox_max,acc,dist_rel,nproxy);
    
    ProxyResultsGamma2(i,2) = toc(start_comp);
    ProxyResultsGamma2(i,5) = NODES_P{9,2};
    ProxyResultsGamma2(i,6) = NODES_P{9,3};
    
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour 3 -- star3 ---

flag_geom='star';
a=0.5; w=10; c=0.6;
params=[a w c];

% BRUTE FORCE COMPRESSION
%----------------------------

disp('Gamma 3: ');
disp('----------');
disp('Brute Force Compression: ');
disp(' ');

for i=1:length(N_vec_BF);
    
    N=N_vec_BF(i);
    nbox_max=nbox_max_vec_BF(i);
    disp(N);
    
    % constructing curve
    [C3] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C3(1,:), C3(2,:), C3(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_BF,Atilde_BF] = OMNICONT_compress_HSS_dsym_brute(C3,nbox_max,acc);
    
    BruteResultsGamma3(i,2) = toc(start_comp);
    BruteResultsGamma3(i,5) = NODES_BF{9,2};
    BruteResultsGamma3(i,6) = NODES_BF{9,3};
    
    btilde=zeros(size(Atilde_BF,1),1);
    btilde(1:N)=b;
    
    %----umfpack--------
    start_lu=tic;
    [L, U, P, Q, R]=lu(Atilde_BF);
    BruteResultsGamma3(i,3) = toc(start_lu);
    
    start_solve=tic;
    sigma_tilde_BF = Q*(U\(L\(P*(R\btilde))));
    BruteResultsGamma3(i,4) = toc(start_solve);
    sigma_BF=sigma_tilde_BF(1:N);
    
    % Evaluating solution
    % test contour
    [theta,phi]=meshgrid(4*pi/5,linspace(0,2*pi,200));
    x=cos(phi).*sin(theta);
    y=sin(phi).*sin(theta);
    z=cos(theta);
    
    h=2*pi/N;
    
    u_BF=0;
    for j=1:N
        u_BF=u_BF+(h*sigma_BF(j).*C3(19,j)/(2*pi)).*(((x-C3(1,j)).*(C3(7,j))+(y-C3(2,j)).*(C3(8,j))+(z-C3(3,j)).*(C3(9,j))))...
            ./((x-C3(1,j)).^2+(y-C3(2,j)).^2+(z-C3(3,j)).^2);
    end
    
    BruteResultsGamma3(i,1) = norm(abs(u_BF-exact(x,y,z)),inf);
end

% PROXY COMPRESSION
%----------------------------

disp('Proxy Compression: ');
disp(' ');

for i=1:length(N_vec_P);
    N=N_vec_P(i);
    nbox_max=nbox_max_vec_P(i);
    disp(N);
    
    % constructing curve
    [C3] = OMNICONT_construct_contour(flag_geom,params,contour_c, N);
    
    % RHS
    g=exact(C3(1,:), C3(2,:), C3(3,:));
    b=(2*g)';
    
    start_comp=tic;
    [NODES_P,Atilde_P] = OMNICONT_compress_HSS_dsym_green(C3,nbox_max,acc,dist_rel,nproxy);
    
    ProxyResultsGamma3(i,2) = toc(start_comp);
    ProxyResultsGamma3(i,5) = NODES_P{9,2};
    ProxyResultsGamma3(i,6) = NODES_P{9,3};
    
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Results:

% COLORS
b1=[0 0.4470 0.7410];
r1=[0.85 0.325 0.0980];
g1=[0.4660 0.6740 0.1880];

% PLOT CONTOURS
% -------------------------------------------------------------------------
figure;
[theta_a,phi_a]=meshgrid(linspace(0, pi,100),linspace(0,2*pi,200));
xg=cos(phi_a).*sin(theta_a);
yg=sin(phi_a).*sin(theta_a);
zg=cos(theta_a);

surf(xg,yg,zg); shading interp;
axis equal; alpha(0.9);
axis off;
cmap=gray(100);
colormap(cmap(60:100,:));
hold on;
h7=plot3(C1(1,:),C1(2,:),C1(3,:),'Color',b1,'LineWidth',2.5);
h8=plot3(C2(1,:),C2(2,:),C2(3,:),'Color',r1,'LineWidth',2.5);
h9=plot3(C3(1,:),C3(2,:),C3(3,:),'Color',g1,'LineWidth',2.5);
set(gcf, 'Color', 'w');
leg=legend([h7 h8 h9], '\Gamma_1','\Gamma_2','\Gamma_3','Orientation','vertical');
set(leg, 'FontSize',14);
set(leg, 'Position',[0.82 0.4 0.17 0.199]);
view([-1.1 85.2]);

% **** for exporting figures:
% addpath(export_fig-master)
% export_fig Ex4-1Domains -m4 -png

% ERROR PLOTS
% -------------------------------------------------------------------------
figure; hold on
%C1 brute
h1_BF = loglog(N_vec_BF, BruteResultsGamma1(:,1), '--','Color',b1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma1(:,1), '.' ,'Color','k','MarkerSize',9);
%C1 proxy
h1_P = loglog(N_vec_P, ProxyResultsGamma1(:,1),'Color',b1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma1(:,1), '.','Color','k', 'MarkerSize',9);
%C2 brute
h2_BF = loglog(N_vec_BF, BruteResultsGamma2(:,1), '--','Color',r1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma2(:,1), '.' ,'Color','k','MarkerSize',9);
%C2 proxy
h2_P = loglog(N_vec_P, ProxyResultsGamma2(:,1),'Color',r1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma2(:,1), '.','Color','k', 'MarkerSize',9);
%C3 brute
h3_BF = loglog(N_vec_BF, BruteResultsGamma3(:,1), '--','Color',g1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma3(:,1), '.' ,'Color','k','MarkerSize',9);
%C3 proxy
h3_P = loglog(N_vec_P, ProxyResultsGamma3(:,1),'Color',g1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma3(:,1), '.', 'Color','k','MarkerSize',9);

set(gca, 'XScale', 'log', 'YScale', 'log')
title('Error');
ylabel('{||\cdot||}_{\infty}');
xlabel('N');
leg=legend([h1_P h1_BF h2_P h2_BF h3_P h3_BF], '\Gamma_1 proxy', '\Gamma_1 brute', '\Gamma_2 proxy', '\Gamma_2 brute', '\Gamma_3 proxy', '\Gamma_3 brute');
set(gcf, 'Color', 'w');
set(gca, 'FontSize',14);

% **** for exporting figures:
% export_fig Ex4-1Error -pdf

% TIMING PLOTS
% -------------------------------------------------------------------------

% Brute Compression Time
%------------------------
figure; hold on
% C1 brute
loglog(N_vec_BF, BruteResultsGamma1(:,2),'--','Color',b1 ,'LineWidth',2.5);
loglog(N_vec_BF, BruteResultsGamma1(:,2),'.' ,'Color','k','MarkerSize',9);
% C2 brute
loglog(N_vec_BF, BruteResultsGamma2(:,2),'--','Color',r1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma2(:,2),'.' ,'Color','k','MarkerSize',9);
% C3 brute
loglog(N_vec_BF, BruteResultsGamma3(:,2),'--','Color',g1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma3(:,2),'.' ,'Color','k','MarkerSize',9);

% plotting O(N^2) Line:
N_vec_BF_Line=2.^(6:16);
loglog(N_vec_BF_Line, 0.000003*N_vec_BF_Line.^2,'k:','LineWidth',2);
text(2^16, 9000, 'O(N^2)', 'FontSize',14);

set(gca, 'XScale', 'log', 'YScale', 'log')
title('Time for Brute Compression');
ylabel('Time (s)');
xlabel('N');
set(gca, 'FontSize',14);
set(gcf, 'Color', 'w');

% **** for exporting figures:
% export_fig Ex4-1TimeBruteComp -pdf

% Proxy Compression Time
%------------------------
figure; hold on
% C1 proxy
loglog(N_vec_P, ProxyResultsGamma1(:,2),'Color',b1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma1(:,2),'.','Color','k','MarkerSize',9);
% C2 proxy
loglog(N_vec_P, ProxyResultsGamma2(:,2),'Color',r1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma2(:,2),'.','Color','k','MarkerSize',9);
% C3 proxy
loglog(N_vec_P, ProxyResultsGamma3(:,2),'Color',g1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma3(:,2),'.','Color','k','MarkerSize',9);

% plotting O(N) Line:
N_vec_P_Line = 2.^(6:18);
loglog(N_vec_P_Line, 0.0002*N_vec_P_Line,'k:','LineWidth',2);
text(2^18, 40, 'O(N)','FontSize',14);

set(gca, 'XScale', 'log', 'YScale', 'log')
title('Time for Proxy Compression');
ylabel('Time (s)');
xlabel('N');
set(gca, 'FontSize',14);
set(gcf, 'Color','w');

% **** for exporting figures:
% export_fig Ex4-1TimeProxyComp -pdf

% Brute/Proxy LU Decomposition Time
%-----------------------------------
figure; hold on
%C1 brute
h4_BF = loglog(N_vec_BF, BruteResultsGamma1(:,3), '--','Color',b1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma1(:,3), '.','Color' ,'k','MarkerSize',9);
%C1 proxy
h4_P = loglog(N_vec_P, ProxyResultsGamma1(:,3),'Color',b1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma1(:,3), '.','Color','k','MarkerSize',9);
%C2 brute
h5_BF = loglog(N_vec_BF, BruteResultsGamma2(:,3), '--','Color',r1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma2(:,3), '.','Color' ,'k','MarkerSize',9);
%C2 proxy
h5_P = loglog(N_vec_P, ProxyResultsGamma2(:,3), 'Color',r1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma2(:,3), '.','Color','k','MarkerSize',9);
%C3 brute
h6_BF = loglog(N_vec_BF, BruteResultsGamma3(:,3), '--','Color',g1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma3(:,3), '.','Color' ,'k','MarkerSize',9);
%C3 proxy
h6_P = loglog(N_vec_P, ProxyResultsGamma3(:,3), 'Color',g1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma3(:,3), '.','Color','k','MarkerSize',9);

% plotting O(N) Line:
N_vec_Line=2.^(5:18);
loglog(N_vec_Line, 0.0001*N_vec_Line,'k:','LineWidth',2);
text(2^18, 30, 'O(N)','FontSize',14);

set(gca, 'XScale', 'log', 'YScale', 'log');
title('Time for LU');
ylabel('Time (s)');
xlabel('N');
set(gca, 'FontSize',14);
leg2=legend([h4_P h4_BF h5_P h5_BF h6_P h6_BF], '\Gamma_1 proxy', '\Gamma_1 brute', '\Gamma_2 proxy', '\Gamma_2 brute', '\Gamma_3 proxy', '\Gamma_3 brute');
set(leg2, 'Location', 'Best');
set(gcf, 'Color', 'w');

% **** for exporting figures:
% export_fig Ex4-1TimeBruteProxyLU -pdf

% Brute/Proxy Solution Time
%---------------------------
figure; hold on
% C1 brute
loglog(N_vec_BF, BruteResultsGamma1(:,4), '--','Color',b1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma1(:,4), '.' ,'Color','k','MarkerSize',9);
% C1 proxy
loglog(N_vec_P, ProxyResultsGamma1(:,4), 'Color',b1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma1(:,4), '.','Color','k','MarkerSize',9);
% C2 brute
loglog(N_vec_BF, BruteResultsGamma2(:,4), '--','Color',r1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma2(:,4), '.','Color' ,'k','MarkerSize',9);
% C2 proxy
loglog(N_vec_P, ProxyResultsGamma2(:,4), 'Color',r1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma2(:,4), '.','Color','k','MarkerSize',9);
% C3 brute
loglog(N_vec_BF, BruteResultsGamma3(:,4), '--','Color',g1 ,'LineWidth' ,2.5);
loglog(N_vec_BF, BruteResultsGamma3(:,4), '.','Color' ,'k','MarkerSize',9);
% C3 proxy
loglog(N_vec_P, ProxyResultsGamma3(:,4), 'Color',g1,'LineWidth',2);
loglog(N_vec_P, ProxyResultsGamma3(:,4), '.','Color','k','MarkerSize',9);

% plotting O(N) Line:
N_vec_Line=2.^(5:19);
loglog(N_vec_Line, 0.000001*N_vec_Line,'k:','LineWidth',2);
text(2^19, 0.5, 'O(N)','FontSize',14);

set(gca, 'XScale', 'log', 'YScale', 'log');
title('Time for Solve');
ylabel('Time (s)');
xlabel('N');
set(gcf, 'Color', 'w');
set(gca,'Color','w')
set(gca, 'FontSize',14);

% **** for exporting figures:
% export_fig Ex4-1TimeBruteProxySolve -pdf

