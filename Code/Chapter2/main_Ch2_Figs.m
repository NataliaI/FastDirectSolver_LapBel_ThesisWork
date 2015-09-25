% Chapter 2 : Outputs figures 2.1 and 2.2 which show an example of a simply
%             connected domain on the sphere and in the stereographic plane
%--------------------------------------------------------------------------

close all; clear; clc;

% SETTINGS: 
%--------------------------------------------------------------------------

% Boundary: 
flag_geom = 'ellipse'; 
a0 = 0.6; b0 = 0.4;
params = [a0 b0]; 

% Location of contour: 
contour_c = [0 0 1]; 

% number of points:
N = 2^10; 
h = 2*pi/N; 

% Direct Solver Settings:
nbox_max = 2^6;
dist_rel = 1.5; 
nproxy   = 50; 
acc      = 1e-12;

% View for plots: 
az = -23.1; el = 45.2;
v  = [az el];

% colors: 
green3 = [0 0.7 0.2];
blue7 = [0 0.5 0.9];
gray1 = [0.827451 0.827451 0.827451];

%--------------------------------------------------------------------------

% constructing contour:
[C] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);

% grid over sphere for plotting: 
n = 1000;
[th,ph] = meshgrid(linspace(0,pi,n),linspace(0,2*pi,2*n)); 
x = cos(ph).*sin(th);
y = sin(ph).*sin(th); 
z = cos(th);

% evaluating integral to determine which points lie in domain Omega
GridMat=zeros(size(x));
for i = 1:N
    GridMat = GridMat + (h*C(19,i)./(2*pi)).*((x-C(1,i)).*(C(7,i))+(y-C(2,i)).*(C(8,i))+(z-C(3,i)).*(C(9,i)))...
                                           ./((x-C(1,i)).^2+(y-C(2,i)).^2+(z-C(3,i)).^2);
end

GridMat(GridMat >= 0) = 1;
GridMat(GridMat < 0)  = 0; 
        
% Plotting Domain on Sphere (Figure 2.1): 
%--------------------------------------------------------------------------
figure; 
hfig = surf(x,y,z, GridMat);
colormap([green3; blue7]); shading interp;
axis equal; hold on 
camlight
hfig.FaceLighting = 'gouraud';
hfig.AmbientStrength = 0.5;
hfig.DiffuseStrength = 0.7;
hfig.SpecularStrength = 0.2;
hfig.SpecularExponent = 25;
hfig.BackFaceLighting = 'unlit';

% Plotting lines of latitude
[th_lat,ph_lat] = meshgrid(linspace(0,pi,6),linspace(0,2*pi,50)); 
x_lat = cos(ph_lat).*sin(th_lat);
y_lat = sin(ph_lat).*sin(th_lat); 
z_lat = cos(th_lat);

for i = 1:6
    plot3(x_lat(:,i),y_lat(:,i),z_lat(:,i), ':','Color', gray1);
end

% Plotting lines of longitude
[th_long,ph_long] = meshgrid(linspace(0, pi,25),linspace(0,2*pi,12));
x_long = cos(ph_long).*sin(th_long);
y_long = sin(ph_long).*sin(th_long); 
z_long = cos(th_long);

for i = 1:12
    plot3(x_long(i,:),y_long(i,:),z_long(i,:), ':','Color', gray1);
end
view(v);
plot3(C(1,:),C(2,:),C(3,:), 'k', 'LineWidth',2);
axis equal; axis off;
set(gcf, 'Color', 'w');

% *** for exporting
% addpath('export_fig-master'); 
% export_fig 'SimplyConnectedDomain' -png -m4
    
% Plotting Domain in Stereographic Plane (Figure 2.2): 
%--------------------------------------------------------------------------    
figure;        
xi   = (x + 1i*y)./(1 - z); 
xi_C = (C(1,:) + 1i*C(2,:))./(1 - C(3,:));
    
plot(real(xi(GridMat==1)), imag(xi(GridMat==1)),'.','MarkerSize',10,'Color',blue7);
hold on 
plot(real(xi(GridMat==0)),imag(xi(GridMat==0)),'.','MarkerSize',10,'Color',green3);
plot(real(xi_C), imag(xi_C), 'k','LineWidth',2);
axis equal; axis([-5 5 -5 5]);
set(gca, 'FontSize',14);

% black outline for plotting grid
x = [-5:0.1:5];
y = -5*ones(size(x));
plot(x,y,'k','LineWidth',2);

y = 5*ones(size(x));
plot(x,y,'k','LineWidth',2);
    
y = [-5:0.1:5];
x = -5*ones(size(y));
plot(x,y,'k','LineWidth',2);
    
x = 5*ones(size(y));
plot(x,y,'k','LineWidth',2);
    
xlabel('Re(\xi)');
ylabel('Im(\xi)');
set(gcf, 'Color', 'w'); 
    
% plotting lines of latitude and longitude   
for j=1:6
    xi_lat = (x_lat(:,j) + 1i*y_lat(:,j))./(1 - z_lat(:,j));
    plot(real(xi_lat), imag(xi_lat), ':','Color',gray1);
end

for j=1:12
    xi_long = (x_long(j,:) + 1i*y_long(j,:))./(1-z_long(j,:));
    plot(real(xi_long), imag(xi_long),':','Color',gray1);
end
    
% export_fig 'SimplyConnectedStereo' -png -m4
  
    
    
