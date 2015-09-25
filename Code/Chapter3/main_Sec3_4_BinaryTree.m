% Sec 3.4 - example of how a boundary contour is divided based on binary
%           tree
%         - ouputs Figure 3.3
%--------------------------------------------------------------------------

clear; clc; close all;

% SETTINGS
%--------------------------------------------------------------------------

% Boundary: 
flag_geom = 'star'; 
a = 0.3; w = 5; c = 0.6; 
params = [a w c]; 

% Location of contour: 
contour_c_ph = 0; 
contour_c_th = pi/4; 
contour_c = [cos(contour_c_ph)*sin(contour_c_th) ...
             sin(contour_c_ph)*sin(contour_c_th) ...
             cos(contour_c_th)];
         
% number of points 
N = 400; 
h = 2*pi/N;

% View for plots: 
az = 92.8; el = 42.8; v = [az el];
% colors:
blue1=[ 0.117647 0.564706 1];
%--------------------------------------------------------------------------

% Constructing contour: 
[C] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);

% Plotting each level of binary tree: 
% Using fixed values for this example: 

% sphere: 
[th,ph] = meshgrid(linspace(0, pi,100),linspace(0,2*pi,200)); 
x = cos(ph).*sin(th);
y = sin(ph).*sin(th); 
z = cos(th);

% Figure 1 : root of tree 
figure; hold on 
surf(x,y,z);
shading interp;
axis equal; alpha(1);
cmap = gray(100);
colormap(cmap(50:100,:)); 
grid off; axis off; 
set(gcf, 'Color', 'w');

plot3(C(1,:),C(2,:),C(3,:),'k-','LineWidth',2.5); 
view(v);

% *** for exporting:
% addpath('export_fig-master'); 
% export_fig 'ContourRoot' -png -m4

% Figure 2 : coarsest level (level 3 in this example) 
figure; hold on 
surf(x,y,z);
shading interp;
axis equal; alpha(1);
colormap(cmap(50:100,:));
grid off; axis off; 
set(gcf, 'Color', 'w'); 

plot3(C(1,1:200),C(2,1:200),C(3,1:200),'Color', blue1,'LineWidth',2.5); 
plot3(C(1,201:end),C(2,201:end),C(3,201:end),'k-','LineWidth',2.5);
set(gcf, 'Color', 'w');
view(v);

% export_fig 'ContourLevel3' -png -m4


% Figure 3 : level 2 
figure; hold on 
surf(x,y,z);
shading interp;
axis equal; alpha(1);
colormap(cmap(50:100,:));
grid off; axis off; hold on 
set(gcf, 'Color', 'w');

plot3(C(1,1:100),C(2,1:100),C(3,1:100),'Color', blue1,'LineWidth',2.5); 
plot3(C(1,101:200),C(2,101:200),C(3,101:200),'r-','LineWidth',2.5);
plot3(C(1,201:300),C(2,201:300),C(3,201:300),'Color', blue1,'LineWidth',2.5);
plot3(C(1,301:400),C(2,301:400),C(3,301:400),'k-','LineWidth',2.5);
view(v);

% export_fig 'ContourLevel2' -png -m4


% Figure 4: level 1 (finest level) 
figure; hold on
surf(x,y,z);
shading interp;
axis equal; alpha(1);
colormap(cmap(50:100,:));
grid off; axis off; hold on 
set(gcf, 'Color', 'w');

plot3(C(1,1:50),C(2,1:50),C(3,1:50),'Color', blue1,'LineWidth',2.5); 
plot3(C(1,51:100),C(2,51:100),C(3,51:100),'r-','LineWidth',2.5);
plot3(C(1,101:150),C(2,101:150),C(3,101:150),'k-','LineWidth',2.5);
plot3(C(1,151:200),C(2,151:200),C(3,151:200),'Color', blue1,'LineWidth',2.5);
plot3(C(1,201:250),C(2,201:250),C(3,201:250),'r-','LineWidth',2.5); 
plot3(C(1,251:300),C(2,251:300),C(3,251:300),'k-','LineWidth',2.5);
plot3(C(1,301:350),C(2,301:350),C(3,301:350),'Color', blue1,'LineWidth',2.5);
plot3(C(1,351:400),C(2,351:400),C(3,351:400),'r-','LineWidth',2.5);
view(v);

% export_fig 'ContourLevel1' -png -m4

