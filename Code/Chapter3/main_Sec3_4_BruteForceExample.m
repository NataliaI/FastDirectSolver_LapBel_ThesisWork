% Section 3.4 - Brute Force Compression Example 
%-----------------------------------------------
% -This script implements an example of brute force compression 
% -Outputs all of the figures in Section 3.4 showing the skeletons from
%  brute force compression 
%--------------------------------------------------------------------------
clear; clc; close all;

% SETTINGS
%--------------------------------------------------------------------------

% Boundary : 
flag_geom = 'star'; 
a = 0.3; w = 5; c = 0.6;
params = [a w c];

% Location of contour : 
contour_c_ph = 0; 
contour_c_th = pi/4; 
contour_c = [cos(contour_c_ph)*sin(contour_c_th) ...
             sin(contour_c_ph)*sin(contour_c_th) ...
             cos(contour_c_th)];

% number of points 
N = 400; 
h = 2*pi/N;

% Direct Solver Settings: 
nbox_max = 50; % number of points in one block on finest level of recursive 
               % algorithm
acc      = 1e-5; % accuracy of ID 

% View for plots: 
az = 92.8; el = 42.8; v = [az el];

%--------------------------------------------------------------------------

% Constructing contour: 
[C] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);
    
% Plotting surface of sphere:    
[th, ph] = meshgrid(linspace(0, pi,100),linspace(0,2*pi,200)); 
x = cos(ph).*sin(th);
y = sin(ph).*sin(th); 
z = cos(th);
    
figure(1);
surf(x,y,z); shading interp; alpha(0.8); 
hold on 
cmap = gray(100);
colormap(cmap(50:100,:));
axis equal; axis off; 
set(gcf, 'Color', 'w');
view(v);
    
% BRUTE FORCE COMPRESSION 
%--------------------------
[Atilde, ~] = OMNICONT_compress_HSS_dsym_brute_Sec3_4(C,nbox_max,acc,v);
    
    
