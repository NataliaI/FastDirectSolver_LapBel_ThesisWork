% Sec 3.1 - Example of rank-deficient off-diagonal blocks
%         - plots off-diagonal singular values of a 4x4 block matrix
%         - outputs Figure 3.2
%--------------------------------------------------------------------------

clear; clc; close all;

% SETTINGS:
%--------------------------------------------------------------------------

% Boundary:
flag_geom = 'star';
a = 0.3; w = 5; c = 0.6;
params = [a w c];

% number of points
N = 400;
h = 2*pi/N;
% block size:
block_size = 100;

% Location of contour:
contour_c_ph = 0;
contour_c_th = pi/4;
contour_c = [cos(contour_c_ph)*sin(contour_c_th) ...
             sin(contour_c_ph)*sin(contour_c_th) ...
             cos(contour_c_th)];

% View for plots:
az = 92.8; el = 42.8; v = [az el];

%--------------------------------------------------------------------------

% constructing contour:
[C] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);

% exact solution and RHS
% location of singularity:
th_0 = pi/4; ph_0 = pi/20;
x0 = [cos(ph_0)*sin(th_0) ...
      sin(ph_0)*sin(th_0) ...
      cos(th_0)];

exact = @(x,y,z) real(1./(((x+1i*y)./(1-z))-((x0(1)+1i*x0(2))./(1-x0(3)))));
g     = exact(C(1,:), C(2,:), C(3,:));
b     = (2*g)';

% constructing system matrix:
A = OMNICONT_construct_A_diag(C,1:N);

% checking accuracy of solution (just to verify everything is correct)
sigma = A\b;

% Evaluating solution
% test contour
[th,ph] = meshgrid(4*pi/5,linspace(0,2*pi,200));
x = cos(ph).*sin(th);
y = sin(ph).*sin(th);
z = cos(th);

u = 0;
for j = 1:N
    u = u + (h*sigma(j).*C(19,j)/(2*pi)).*(((x-C(1,j)).*(C(7,j))+(y-C(2,j)).*(C(8,j))+(z-C(3,j)).*(C(9,j))))./((x-C(1,j)).^2+(y-C(2,j)).^2+(z-C(3,j)).^2);
end

err = norm(abs(u-exact(x,y,z)),inf);
disp(['max solution error: ' num2str(err)]);

% PLOTTING OFF-DIAGONAL BLOCKS:

Aplot = zeros(N,N);
Aplot(1:block_size, block_size + 1:end) = 1;
Aplot(block_size+1:2*block_size, [1:block_size 2*block_size+1:4*block_size]) = 1;
Aplot(2*block_size+1:3*block_size, [1:2*block_size 3*block_size + 1:4*block_size]) = 1;
Aplot(3*block_size:end, 1:3*block_size) = 1;

spy(Aplot);

% PLOTTING SINGULAR VALUES OF OFF-DIAGONAL BLOCKS:
SV1 = svd(A(1:block_size, block_size+1:end));
SV2 = svd(A(block_size+1:2*block_size, [1:block_size 2*block_size + 1:4*block_size]));
SV3 = svd(A(2*block_size + 1:3*block_size, [1:2*block_size 3*block_size + 1:4*block_size]));
SV4 = svd(A(3*block_size +1 :end, 1:3*block_size));

figure; hold on
semilogy(SV1,'k.');
semilogy(SV2,'k.');
semilogy(SV3,'k.');
semilogy(SV4,'k.');
set(gca, 'XScale', 'linear', 'YScale', 'log')
set(gcf, 'Color', 'w');
set(gca, 'FontSize',14);
xlabel('Ordered Index');
ylabel('Singular Values');

% ** for exporting: 
% addpath('export_fig-master'); 
% export_fig 'OffDiagSingValues' -pdf 

% PLOTTING BOUNDARY:
[th,ph] = meshgrid(linspace(0, pi,100),linspace(0,2*pi,200));
x = cos(ph).*sin(th);
y = sin(ph).*sin(th);
z = cos(th);

figure; hold on
surf(x,y,z); 
shading interp; alpha(1); 
cmap = gray(100);
colormap(cmap(50:100,:));

plot3(C(1,:),C(2,:),C(3,:),'k.','MarkerSize',8);
axis equal; axis off;
set(gcf, 'Color', 'w');
view(v);
hold on

% export_fig 'OffDiagExContour' -pdf 

