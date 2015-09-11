function [] = main_Ex6_6SDC_CapConvergence()
%--------------------------------------------------------------------------
% Ex 6-6 : Performance of SDC for obtaining point vortex trajectories
%          in a spherical cap domain
% Generates convergence plots in Figure 6.8 and time plot in Figure 6.9
% written as a function so that helper functions can be included in file
%--------------------------------------------------------------------------

close all; clear; clc;

% Settings:
%--------------------------------------------------------------------------

% Location of Cap and Point Vortex on Cap:
%----------------------------------------
% phi   - azimuthal [0,2pi]
% theta - elevation [0, pi]

th_c = pi/6;                % spherical cap
th_1 = pi/20; ph_1=pi/20;   % point vortex on cap
x1   = [cos(ph_1)*sin(th_1) sin(ph_1)*sin(th_1) cos(th_1)];

% Boundary:
%----------
gamma = 1;        % vortex strength

flag_geom = 'cap';
params    = [sin(th_c) th_c];
contour_c = [0 0 1];

N = 100;
disp(['N = ' num2str(N)]);
nbox_max = 10;
h = 2*pi/N;

% Direct Solver Settings:
%------------------------
acc      = 1e-12;
dist_rel = 1.5;
nproxy   = 50;

% Settings for Plots:
%--------------------
% colors:
b1=[0 0.4470 0.7410];
r1=[0.85 0.325 0.0980];
%--------------------------------------------------------------------------

% constructing curve
[C] = OMNICONT_construct_contour(flag_geom, params, contour_c, N);

% boundary data
boundary = @(x,y,z,x0,x1) (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2))...
    -(gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2));

% COMPRESSING Atilde
%---------------------
time_comp=tic;
[NODES, Atilde] = OMNICONT_compress_HSS_dsym_green(C,nbox_max,acc,dist_rel,nproxy);
disp(['Time for compression: ' num2str(toc(time_comp))]);

btilde=zeros(size(Atilde,1),1);

%----umfpack--------------
[L, U, P, Q, R]=lu(Atilde);

%--------------------------------------------------------------------------
% SDC

% Settings:
%------------
J  = 9;
Mp = 10;

th_0 = pi/2;  ph_0 = 0;      % initial condition (vortex position)
ICx0 = [cos(ph_0)*sin(th_0) sin(ph_0)*sin(th_0) cos(th_0)];

% Obtaining initial velocity (for exact solution)
sigma = LOCAL_get_sigma(boundary, C, ICx0, x1, btilde, L, U, P, Q, R);
vel_initial = LOCAL_eval_F(1, ICx0, x1, gamma, sigma, C, h);
norm_v = norm(vel_initial);
T = (2*pi)/norm_v; % time for full trajectory
T = T/2;           % just evaluating for half of trajectory

count=1;

% Convergence rate for varying N
%--------------------------------------------------------------------------
npanel_vec = 1:20;

time_SDC = zeros(length(npanel_vec),1);
errM     = zeros(length(npanel_vec),1);

disp(' ')
disp('Finding convergence rate for varying N: ');
disp(' ');
disp(['Number of points on each panel = ' num2str(Mp)]);
disp(['Number of correction steps = ' num2str(J)]);
disp(' ');

for npanels = npanel_vec
    disp(['Number of panels ' num2str(npanels)]);
    
    start_SDC=tic;
    [sol, time, ~] = SDC(ICx0, x1, gamma, C, T, J, Mp, npanels, N, boundary, btilde, L, U, P, Q, R);
    time_SDC(count) = toc(start_SDC);
    
    % exact solution
    xexact=cos(norm_v*time);
    yexact=-sin(norm_v*time);
    zexact=zeros(length(time),1);
    
    errM(count)=norm([xexact yexact zexact]-sol,inf);
    count=count+1;
end

M = Mp*npanel_vec;

% Plotting Convergence Rate
figure; hold on
loglog(M, errM,'Color',b1,'LineWidth',2);
loglog(M, errM,'k.','MarkerSize',9);
loglog(M, 4200*(M).^(-(J+1)),'k--', 'LineWidth',2);
loglog(M, 50*(M).^(-(J+(1/2))),'--','Color',r1,'LineWidth',2);

set(gca, 'XScale', 'log', 'YScale', 'log')
text(100,10^(-16), 'O(M^{-(J+1)})','FontSize',14);
text(40,10^(-17), 'O(M^{-(J+1/2)})','FontSize',14);
title(['Convergence rate for varying M and fixed J = ' num2str(J)]);
xlabel('M');
ylabel('||\cdot||_{\infty}');
set(gca, 'FontSize',14);
set(gcf, 'Color', 'w');
%*** for exporting: 
% addpath('export_fig-master'); 
% export_fig Ex6-6CapSDCConvergenceVaryMfixJ -pdf 

% Plotting Solution
[th, ph]=meshgrid(linspace(0,pi,100), linspace(0,2*pi,200));
x = cos(ph).*sin(th);
y = sin(ph).*sin(th);
z = cos(th);

figure; hold on
surf(x,y,z); shading interp;
h1 = plot3(sol(:,1),sol(:,2), sol(:,3),'b', 'LineWidth',2);
h2 = plot3(C(1,:),C(2,:),C(3,:), 'k', 'LineWidth',2);
axis equal; axis off;
set(gcf, 'Color','w');
leg = legend([h1 h2], 'solution', 'boundary');
set(leg, 'FontSize',14);
az=-23.5;  el=32; v=[az el];
view(v);

% Plotting Timings vs M
figure;
hold on
loglog(M, time_SDC, 'Color', b1, 'LineWidth',2);
loglog(M, time_SDC,'k.', 'MarkerSize',9);
loglog(M, 0.02*M,'k--');

text(20,0.8,'O(M)','FontSize',16);
set(gca, 'XScale', 'log', 'YScale', 'log')
ylabel('Time (s)');
xlabel('M');
set(gca, 'FontSize',14);
set(gcf, 'Color', 'w');
% ***
% export_fig Ex6-6CapSDCTime -pdf 

% Convergence rate for varying J
%--------------------------------------------------------------------------

npanels = 10;

disp(' ')
disp('Finding convergence rate for varying J: ');
disp(' ');
disp(['Number of points on each panel = ' num2str(Mp)]);
disp(['Number of panels (fixed) = ' num2str(npanels)]);
disp(['Number of correction steps = ' num2str(J)]);
disp(' ');

M = Mp*npanels;
pwidth = T/npanels;     % panel width
a = 0; Tp = a+pwidth;   % panel interval [a, Tp]
sol = zeros(M,3);       % for storing solutions for each panel
t  = zeros(M,1);
x0 = ICx0;

[tp]   = chebpts(Mp);      % grid for one panel on [-1,1]

dt=zeros(length(tp)-1,1);  % step sizes
for i=1:Mp-1
    dt(i)=tp(i+1)-tp(i);
end

for p=1:npanels
    x      = zeros(Mp,3);
    
    % Forward Euler (working on t in [-1,1])
    x(1,:) = ICx0;              % storing solution (values of x0(t))
    
    sigma_mat = zeros(N,Mp); % stores solutions for sigma for each time step
    for step = 1:Mp-1
        sigma = LOCAL_get_sigma(boundary, C, x0, x1, btilde, L, U, P, Q, R);
        sigma_mat(:,step) = sigma;
        
        vel = LOCAL_eval_F(((Tp-a)/2)*tp(step)+(Tp+a)/2, x0, x1, gamma, sigma, C, h); % velocity on original interval [a,T]
        
        x(step+1,:)      = x(step,:)+dt(step)*((Tp-a)/2)*vel;   % stretches step h to interval [a,T]
        
        x0 = x(step+1,:);  % on original interval
    end
    
    sigma = LOCAL_get_sigma(boundary, C, x0, x1, btilde, L, U, P, Q, R);  % plug into matrix
    sigma_mat(:,step+1) = sigma;
    
    S = cumsummat(Mp);  % spectral integration matrix
    
    errJ=zeros(J,1);
    
    for j = 1:J   % correction steps
        
        disp(['Correction step j = ' num2str(j)]);
        
        % compute residual function (for each t)
        vel = LOCAL_eval_F_mat((Tp/2)*(1+tp), x, x1, gamma, sigma_mat, C, h);
        residual = S*((Tp-a)/2)*vel-x+[ICx0(1)*ones(Mp,1) ICx0(2)*ones(Mp,1) ICx0(3)*ones(Mp,1)];
        
        % Forward Euler to compute error (delta)
        %---------------------------------------
        delta = zeros(Mp,3);  % error
        delta(1,:) = residual(1,:);
        for step = 1:Mp-1
            
            % computing G(t(step), delta(step)) :
            
            % F(t(step), x(step)+delta(step))
            sigmaF1 = LOCAL_get_sigma(boundary, C, x(step,:) + delta(step,:), x1, btilde, L, U, P, Q, R);
            F1 = LOCAL_eval_F(((Tp-a)/2)*tp(step)+(Tp+a)/2, x(step,:) + delta(step,:), x1, gamma, sigmaF1, C, h);
            
            % F(t(step), x(step))
            sigmaF2 = LOCAL_get_sigma(boundary, C, x(step,:), x1, btilde, L, U, P, Q, R);
            F2 = LOCAL_eval_F(((Tp-a)/2)*tp(step)+(Tp+a)/2, x(step,:), x1, gamma, sigmaF2, C, h);
            
            % delta(step+1) = delta(step) + dt(step)*G(t(step),delta(step))
            %                 + residual(step+1) - residual(step)
            delta(step+1,:) = delta(step,:) + dt(step)*(((Tp-a)/2)*F1 - ((Tp-a)/2)*F2)...
                + residual(step+1,:) - residual(step,:);
        end
        
        x=x+delta; % updating solution
        
        % recomputing density sigma
        sigma_mat=zeros(N, length(x));
        for k=1:length(x)
            sigma_mat(:,k) = LOCAL_get_sigma(boundary, C, x(k,:), x1, btilde, L, U, P, Q, R);
        end
        
        % exact solution
        xexact = cos(norm_v*(((Tp-a)/2)*tp+(Tp+a)/2));
        yexact = -sin(norm_v*(((Tp-a)/2)*tp+(Tp+a)/2));
        zexact = zeros(length(tp),1);
        
        errJ(j) = norm([xexact yexact zexact]-x,inf); % error for each correction step
    end
    
    % Plotting convergence of error for varying J
    figure(4); hold on
    semilogy(1:J, errJ,'Color',b1,'LineWidth',2);
    semilogy(1:J, errJ,'k.','MarkerSize',9);
    semilogy(1:J, 0.05*exp(-6*(1:J)),'k--','LineWidth',2);
    hold on
    set(gca, 'XScale', 'linear', 'YScale', 'log')
    if p==npanels
        text(3.3, 10^(-9), 'exp(-{(J+1)})','FontSize',14);
    end
    title(['Convergence Rate for Varying J, Fixed M = ',num2str(M)]);
    xlabel('Correction Steps j');
    ylabel('||\cdot||_{\infty}');
    set(gca, 'FontSize',14);
    set(gcf, 'Color', 'w');
    % ****
    % export_fig Ex6-6CapSDCConvergenceVaryJfixM -pdf 
    
    t((p-1)*Mp+1: p*Mp, :) = (Tp-a)/2*tp+(Tp+a)/2;
    sol((p-1)*Mp+1: p*Mp,:) = x;
    
    % initial condition for next panel
    ICx0 = [x(end,1) x(end,2) x(end,3)];
    
    a = Tp;
    Tp = a+pwidth;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sigma = LOCAL_get_sigma(boundary, C, x0, x1, btilde, L, U, P, Q, R)
%--------------------------------------------------------------------------
% This function solves the system Atilde*sigma_tilde = btilde for each new
% vortex position x0
%--------------------------------------------------------------------------
N = size(C,2);
g = boundary(C(1,:),C(2,:), C(3,:), x0,x1);
b = 2*g';
btilde(1:N) = b;
sigma_tilde = Q*(U\(L\(P*(R\btilde))));
sigma = sigma_tilde(1:N);

return

function [vel] = LOCAL_eval_F(t, x0, x1, gamma, sigma, C, h)
%--------------------------------------------------------------------------
% This function evaluates the right hand side of the IVP for each
% vortex position x0(t)
% dx0(t)/dt = F(t,x0(t))
% F(t,x0(t)) = (grad_x \hat{psi}(x,x0) x e_r)|x=x0 = (velocity of vortex)
% F(t,x0(t)) = velocity = (u,v,w)
% Details in Section 5.3
%--------------------------------------------------------------------------
u = sum((h*sigma'.*C(19,:)./(2*pi)).*(C(7,:)./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2)-(C(7,:).*(x0(1)-C(1,:))+...
    C(8,:).*(x0(2)-C(2,:))+C(9,:).*(x0(3)-C(3,:))).*(2*x0(1)-2*C(1,:))./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2).^2))+...
    gamma/(2*pi)*(x0(1)-x1(1))/((x0(1)-x1(1))^2+(x0(2)-x1(2))^2+(x0(3)-x1(3))^2);

v = sum((h*sigma'.*C(19,:)./(2*pi)).*(C(8,:)./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2)-(C(7,:).*(x0(1)-C(1,:))+...
    C(8,:).*(x0(2)-C(2,:))+C(9,:).*(x0(3)-C(3,:))).*(2*x0(2)-2*C(2,:))./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2).^2))+...
    gamma/(2*pi)*(x0(2)-x1(2))/((x0(1)-x1(1))^2+(x0(2)-x1(2))^2+(x0(3)-x1(3))^2);

w = sum((h*sigma'.*C(19,:)./(2*pi)).*(C(9,:)./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2)-(C(7,:).*(x0(1)-C(1,:))+...
    C(8,:).*(x0(2)-C(2,:))+C(9,:).*(x0(3)-C(3,:))).*(2*x0(3)-2*C(3,:))./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2).^2))+...
    gamma/(2*pi)*(x0(3)-x1(3))/((x0(1)-x1(1))^2+(x0(2)-x1(2))^2+(x0(3)-x1(3))^2);

vel = [v*x0(3)-w*x0(2) -(u*x0(3)-w*x0(1)) u*x0(2)-v*x0(1)];
return

function [vel] = LOCAL_eval_F_mat(t,x0_vec, x1, gamma, sigma_mat, C, h)
%--------------------------------------------------------------------------
% This function evaluates the right hand side of the IVP, F(t,x0(t)) for a
% vector of vortex positions x0. It takes a matrix of sigma values as input
% corresponding to each position x0.
%--------------------------------------------------------------------------

vel = zeros(size(x0_vec,1),3);

for k=1:size(x0_vec,1)
    sigma = sigma_mat(:,k);
    x0 = x0_vec(k,:);
    u = sum((h*sigma'.*C(19,:)./(2*pi)).*(C(7,:)./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2)-(C(7,:).*(x0(1)-C(1,:))+...
        C(8,:).*(x0(2)-C(2,:))+C(9,:).*(x0(3)-C(3,:))).*(2*x0(1)-2*C(1,:))./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2).^2))+...
        gamma/(2*pi)*(x0(1)-x1(1))/((x0(1)-x1(1))^2+(x0(2)-x1(2))^2+(x0(3)-x1(3))^2);
    
    v = sum((h*sigma'.*C(19,:)./(2*pi)).*(C(8,:)./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2)-(C(7,:).*(x0(1)-C(1,:))+...
        C(8,:).*(x0(2)-C(2,:))+C(9,:).*(x0(3)-C(3,:))).*(2*x0(2)-2*C(2,:))./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2).^2))+...
        gamma/(2*pi)*(x0(2)-x1(2))/((x0(1)-x1(1))^2+(x0(2)-x1(2))^2+(x0(3)-x1(3))^2);
    
    w = sum((h*sigma'.*C(19,:)./(2*pi)).*(C(9,:)./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2)-(C(7,:).*(x0(1)-C(1,:))+...
        C(8,:).*(x0(2)-C(2,:))+C(9,:).*(x0(3)-C(3,:))).*(2*x0(3)-2*C(3,:))./((x0(1)-C(1,:)).^2+(x0(2)-C(2,:)).^2+(x0(3)-C(3,:)).^2).^2))+...
        gamma/(2*pi)*(x0(3)-x1(3))/((x0(1)-x1(1))^2+(x0(2)-x1(2))^2+(x0(3)-x1(3))^2);
    
    vel(k,:)=[v*x0(3)-w*x0(2) -(u*x0(3)-w*x0(1)) u*x0(2)-v*x0(1)];
end

