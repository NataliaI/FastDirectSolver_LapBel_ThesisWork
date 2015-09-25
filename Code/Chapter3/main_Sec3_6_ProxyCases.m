% Section 3.6 - Proxy Point Compression on the Sphere
%             - examples of different cases on pages 45-46
%             - outputs Figures 3.11 and 3.12
%---------------------------------------------------------

clear; clc; close all;

% SETTINGS
%--------------------------------------------------------------------------

% Boundary: 
flag_geom = 'cap'; 
th_c      = pi/2; 
params    = [sin(pi/2) pi/2]; 

% Location of contour: 
contour_c = [0 0 1]; 

% number of points
N = 400; 
h = 2*pi/N; 

% Direct Solver Settings
nbox_max = 50;   % number of points in one block on finest level of recursive 
                 % algorithm
acc      = 1e-15; % accuracy of ID 
dist_rel = 1.5;  % ratio which determines how proxy circle is scaled
nproxy   = 50; 


% View for plots: 
az = 65.6; el = 10; v = [az el];
%--------------------------------------------------------------------------

% Constructing contour: 
[C] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);

% Plotting surface of sphere:    
[th, ph] = meshgrid(linspace(0, pi,100),linspace(0,2*pi,200)); 
x = cos(ph).*sin(th);
y = sin(ph).*sin(th); 
z = cos(th);

figure(1);
surf(x,y,z); shading interp; alpha(0.5); 
hold on 
cmap = gray(100);
colormap(cmap(50:100,:));
axis equal; axis off; 
set(gcf, 'Color', 'w');
view(v);

% PROXY POINT COMPRESSION 
%--------------------------
[Atilde, ~] = OMNICONT_compress_HSS_dsym_green_Sec3_6Cases(C,nbox_max,acc,dist_rel, nproxy,v);
