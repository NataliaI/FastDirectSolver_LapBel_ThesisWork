% Sec 3.5 - 2D Proxy Points
%         - example of how proxy points work in 2D
%         - plots near field, far field and proxy circle 
%         - outputs Figures 3.7, 3.8 and 3.9
%--------------------------------------------------------------------------

close all; clear; clc; 

% SETTINGS: 
%--------------------------------------------------------------------------

% Boundary: 
flag_geom = 'star'; 
a = 0.5; w = 2; c = 0.2; 
params = [a w c]; 

% number of points 
N = 500; 

% colors for plots: 
green5 = [ 0.133333 0.545098 0.133333];
gray1  = [0.55 0.55 0.55];

%--------------------------------------------------------------------------
% constructing contour: 
t = linspace(0,2*pi*(1-1/N), N); 
R = c*(1+a*cos(w*t)); 

Cx = R.*cos(t); 
Cy = R.*sin(t); 

figure(1); hold on  % Figure 3.8 
                    % plots Gamma_tau, Gamma_near, Gamma_far and
                    % Gamma_proxy
plot(Cx(110:170), Cy(110:170),'r', 'LineWidth', 2.8);
plot(Cx(79:188),  Cy(79:188) , 'k', 'LineWidth', 2.8);
plot(Cx(1:78),    Cy(1:78)   ,'Color', gray1,'LineWidth', 2.8);
plot(Cx(189:end), Cy(189:end),'Color',gray1,'LineWidth', 2.8);
plot(Cx(110:170), Cy(110:170),'r', 'LineWidth',2.8);

text(-0.03, 0.13, '\Gamma_\tau', 'FontSize',18);
text(-0.13, 0.11, '\Gamma_{near}', 'FontSize',18);
text(-0.1, 0.25, '\Gamma_{proxy}', 'FontSize',18);
text(0.2, -0.08, '\Gamma_{far}', 'FontSize',18);

% proxy circle: 
Cpx = 0.12*cos(t) - 0.025;
Cpy = 0.12*sin(t) + 0.11; 

plot(Cpx, Cpy, 'Color', green5, 'LineWidth', 2.8);
axis equal; 
xlimits = get(gca, 'XLim');
ylimits = get(gca, 'YLim');
axis off; 
set(gcf, 'Color', 'w');

% *** for exporting
% addpath('export_fig-master'); 
% export_fig '2DGammaTauNearFarProxy' -pdf -transparent

figure(2); hold on  % Figure 3.7
                    % plots Gamma_tau and Gamma_tau complement
plot(Cx, Cy,'k','LineWidth',2.8);
plot(Cx(110:170), Cy(110:170),'r', 'LineWidth', 2.8);

text(-0.03, 0.13, '\Gamma_\tau', 'FontSize',18);
text(0.2, -0.08, '\Gamma_\tau^c', 'FontSize',18);
xlim(xlimits);
ylim(ylimits);
axis off; 
set(gcf, 'Color', 'w');

% export_fig '2DGammaTauGammaTauC' -pdf -transparent
 
figure(3); hold on  % Figure 3.9
                    % plots Gamma_tau, Gamma_far and Gamma_proxy
plot(Cx(1:78), Cy(1:78),'Color', gray1,'LineWidth',2.8);
plot(Cx(189:end), Cy(189:end),'Color',gray1,'LineWidth',2.8);
plot(Cx(110:170), Cy(110:170),'r', 'LineWidth',2.8);

Cpx = 0.12*cos(t) - 0.025;
Cpy = 0.12*sin(t) + 0.11; 

plot(Cpx, Cpy, 'Color', green5,'LineWidth', 2.8);
axis equal; 
xlimits = get(gca, 'XLim');
ylimits = get(gca, 'YLim');
axis off; set(gcf, 'Color', 'w');

text(-0.03, 0.13, '\Gamma_\tau', 'FontSize',18);
text(-0.1, 0.25, '\Gamma_{proxy}', 'FontSize',18);
text(0.2, -0.08, '\Gamma_{far}', 'FontSize',18);

% export_fig '2DGammaTauFarProxy' -pdf -transparent

