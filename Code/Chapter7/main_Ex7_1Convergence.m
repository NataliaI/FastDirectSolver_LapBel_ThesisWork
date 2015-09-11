%--------------------------------------------------------------------------
% Chapter 7 : A fast direct solver for a multiply connected domain
%           : This script solves the Laplace-Beltrami equation for a domain
%             with two star-shaped islands
%           : nbox_max is set so that the recursive levels always stay on
%             one contour at a time
%           : Produces a convergence plot of the error - Figure 7.1b)
%--------------------------------------------------------------------------

close all; clear; clc;

% Settings for direct solver:
%----------------------------

% number of points/recursive levels :
N_vec        = 2.^(3:8);
nbox_max_vec = [2 2^2 2^3 2^3 2^4 2^4];

err          = zeros(size(N_vec));

acc      = 1e-10;  % accuracy of ID
dist_rel = 1.6;    % ratio which determines how proxy circle is scaled
nproxy   = 60;     % number of points on proxy cap

% Location of all contours
contour_c_ph = [pi/6; pi/8];
contour_c_th = [pi/6; pi/2];
contour_c    = [cos(contour_c_ph).*sin(contour_c_th) ...
                sin(contour_c_ph).*sin(contour_c_th) ...
                cos(contour_c_th)];

% contour parameters
flag_geom = 'star';
params = [0.2 5 0.4;
          0.2 5 0.3];
M = size(params,1);

exact = @(x,y,z) (1/2)*(real(1./(((x+1i*y)./(1-z))-(contour_c(1,1)+1i*contour_c(1,2))./(1-contour_c(1,3)))))...
                +(1/2)*(real(1./(((x+1i*y)./(1-z))-(contour_c(2,1)+1i*contour_c(2,2))./(1-contour_c(2,3)))));

for i = 1:length(N_vec)
    
    N        = N_vec(i);
    nbox_max = nbox_max_vec(i);
    h        = 2*pi/N;
    disp(['N = ' num2str(N)]);
    disp(['nbox_max = ' num2str(nbox_max)]);
    disp(' '); 
    
    % Constructing contour:
    [C,t]  = OMNICONT_construct_contour_MC(flag_geom, params, contour_c, N);
    
    % compressing system
    [Atilde] = OMNICONT_compress_HSS_dsym_green_MC(C,nbox_max,acc,dist_rel, nproxy,N,M,contour_c);
    
    % RHS
    g = exact(C(1,:),C(2,:), C(3,:));
    b = 2*g';
    btilde = zeros(size(Atilde,1),1);
    btilde(1:N*M) = b;
    
    %----umfpack-----------------
    [L, U, P, Q, R] = lu(Atilde);
    sigma_tilde = Q*(U\(L\(P*(R\btilde))));
    sigma = sigma_tilde(1:M*N);
    
    % Evaluating solution
    % (test contour is well away from boundary)
    [theta,phi] = meshgrid(linspace(5*pi/6, 5*pi/6,1),linspace(0,2*pi,100));
    x = cos(phi).*sin(theta);
    y = sin(phi).*sin(theta);
    z = cos(theta);
    
    u = 0;
    for k = 1:M % finding coefficients of log terms : A_k
        A(k) = h*sigma((k-1)*N+1:k*N)'*C(19,(k-1)*N+1:k*N)';
    end
    
    for k = 1:N*M
        u = u+(h*sigma(k).*C(19,k)/(2*pi)).*(((x-C(1,k)).*(C(7,k))+(y-C(2,k)).*(C(8,k))+(z-C(3,k)).*(C(9,k))))./((x-C(1,k)).^2+(y-C(2,k)).^2+(z-C(3,k)).^2);
    end
    for k = 2:M % adding on sources
        u = u + A(k)*((-1/(2*pi))*(log(sqrt((x-contour_c(1,1)).^2+(y-contour_c(1,2)).^2+(z-contour_c(1,3)).^2)/4))...
            +(1/(2*pi))*log(sqrt((x-contour_c(2,1)).^2+(y-contour_c(2,2)).^2+(z-contour_c(2,3)).^2)/4));
    end
    
    err(i) = norm(abs(u - exact(x, y, z)), inf);
    
end

figure;
b1 = [0 0.4470 0.7410];
semilogy(M*N_vec, err,'Color',b1,'LineWidth',2);
hold on
semilogy(M*N_vec, err,'.', 'Color', 'k','MarkerSize',10);
set(gca,'FontSize',14);
xlabel('M*N');
ylabel('{|| \cdot ||}_\infty');
set(gcf, 'Color', 'w');
set(gca,'Color','w') ;
title('Error');

%*** for exporting:
% addpath('export_fig-master');
% export_fig Ex7-1TwoStarsError -pdf




