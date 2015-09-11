%--------------------------------------------------------------------------
% Ex 4-1 Results
% This script loads results generated for Ex4-1 in the thesis and outputs the
% corresponding figures
%--------------------------------------------------------------------------

clear; clc; close all;
load Ex4-1Results.mat

% To export figures: 
% addpath('export_fig-master')
% uncomment export_fig statements at the end of each figure 

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


