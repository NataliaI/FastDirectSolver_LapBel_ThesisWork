%--------------------------------------------------------------------------
% Ex 4-2 Results
% This script loads results generated for Ex4-2 in the thesis and outputs the
% corresponding figures
%--------------------------------------------------------------------------

clear; clc; close all;
load Ex4-2Results.mat

% To export figures: 
% addpath('export_fig-master')
% uncomment export_fig statements at the end of each figure 

% Plotting Results: 
%--------------------------------------------------------------------------

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

h1 = plot3(C1(1,:),C1(2,:),C1(3,:),'Color',r1,'LineWidth',2); 
h2 = plot3(C4(1,:),C4(2,:),C4(3,:),'Color',b1,'LineWidth',2);
set(gcf, 'Color', 'w');
leg=legend([h1 h2], '\Gamma_1','\Gamma_4');
set(leg, 'FontSize',14);
set(leg, 'Position',[0.76 0.4 0.17 0.199]);
view([az el]);

% **** for exporting figures:
% export_fig Ex4-2Gamma1and4 -m4 -png

% Stars: Gamma 2 and Gamma 3
figure; hold on 
surf(xg,yg,zg); shading interp;
axis equal; alpha(0.9);
axis off; 
cmap=gray(100);
colormap(cmap(50:100,:)); 

h3 = plot3(C2(1,:),C2(2,:),C2(3,:),'Color',r1,'LineWidth',2);
h4 = plot3(C3(1,:),C3(2,:),C3(3,:),'Color',b1,'LineWidth',2);

set(gcf, 'Color', 'w');
leg2=legend([h3 h4], '\Gamma_2','\Gamma_3');
set(leg2, 'FontSize',14);
set(leg2, 'Position',[0.76 0.4 0.17 0.199]);
view([az el]);

% **** for exporting figures:
% export_fig Ex4-2Gamma2and3 -m4 -png

% Printing Table Values
%--------------------------------------------------------------------------
disp('Gamma 1 Results: '); 
disp('-----------------');
disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
disp(ProxyResultsGamma1)
disp(' '); 

disp('Gamma 2 Results: '); 
disp('-----------------');
disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
disp(ProxyResultsGamma2)
disp(' '); 

disp('Gamma 3 Results: '); 
disp('-----------------');
disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
disp(ProxyResultsGamma3)
disp(' ');

disp('Gamma 4 Results: '); 
disp('-----------------');
disp(['    Error' '     Time_c' '    Time_lu' '   Time_s' '    k_p' '          K']);
disp(ProxyResultsGamma4)


