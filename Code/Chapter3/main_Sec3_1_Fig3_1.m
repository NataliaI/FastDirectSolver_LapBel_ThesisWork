% Sec 3.1 - Outputs Figure 3.1
%         - divides sources along boundary according to their
%           self-interaction, and the near and far field
%         - only plots one example for a section of the boundary, Gamma_tau
%           on the finest level
%--------------------------------------------------------------------------

function [] = main_Sec3_1_Fig3_1()

clear; clc; close all;

% SETTINGS:
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
N = 800; 

% Direct Solver Settings:
nbox_max = 100 ;  
dist_rel = 1.5;
nproxy   = 50; 

% view for plots:
az = 92.8; el = 42.8; v = [az el];

% colors:
gray1  = [0.55 0.55 0.55];

%--------------------------------------------------------------------------
% constructing contour:
[Ctemp] = OMNICONT_construct_contour(flag_geom,params,contour_c,N);
C = [Ctemp(:,500:end) Ctemp(:,1:499)]; % rearranging indices for nicer plots

% Plotting surface of sphere:
[th, ph] = meshgrid(linspace(0, pi,100),linspace(0,2*pi,200));
x = cos(ph).*sin(th);
y = sin(ph).*sin(th);
z = cos(th);

figure(1); hold on
surf(x,y,z); shading interp;
cmap = gray(100);
colormap(cmap(50:100,:));
plot3(C(1,:),C(2,:),C(3,:),'.','MarkerSize',9,'Color',gray1);
axis equal; axis off;
set(gcf, 'Color', 'w');
alpha(1); view(v);

% Constructing binary/recursive tree
% taking sections from proxy compression code to identify Gamma_tau, near/far-field

NODES  = LOCAL_get_tree(C,nbox_max);
nboxes = size(NODES,2);

% Constructing the list of neighbors of any cell.
LIST_NEI = LOCAL_construct_potential_neighborlist(C,NODES,dist_rel);

ibox = nboxes; % box on finest level
indskel = NODES{6,ibox} - 1 + (1:NODES{7,ibox});

% Determine a spherical cap that circumscribes Gamma_tau (indskel):
[xxc, max_dist] = LOCAL_get_circum_circle(C(:,sort(indskel)));

dist = dist_rel*max_dist; % distance between xxc and a point on spherical cap

[Cproxy] = LOCAL_construct_cap(xxc, dist, nproxy);

% Obtaining neighbours of Gamma_tau :
indskel_maybeinside = [];
for jbox = LIST_NEI{ibox}
    if (NODES{5,jbox}==0)
        indskel_maybeinside = [indskel_maybeinside,NODES{6,jbox}-1+(1:NODES{7,jbox})];
    else
        jbox_son1 = NODES{4,jbox}(1);
        jbox_son2 = NODES{4,jbox}(2);
        indskel_maybeinside = [indskel_maybeinside,NODES{11,jbox_son1},NODES{11,jbox_son2}];
    end
end

z_axis = [xxc(1) xxc(2) xxc(3)]; % creating axis centered at xxc
                                 % (taken as north pole)
                                 % obtaining z-coordinate of proxy cap
zcoord_proxy = Cproxy(1,1)*z_axis(1)+Cproxy(2,1)*z_axis(2)+Cproxy(3,1)*z_axis(3);

% obtaining z-coordinate of an arbitrary point on Gamma_tau
% (used to determine orientation of proxy cap)
zcoord_indskel = C(1,indskel(1))*z_axis(1)+C(2,indskel(1))*z_axis(2)+C(3,indskel(1))*z_axis(3);

if zcoord_indskel>zcoord_proxy
    direction = 'above';  % Gamma_tau above proxy circle
else
    direction = 'below';  % Gamma_tau below proxy circle
end

% obtaining z-coordinates of neighbours of Gamma_tau
zcoord_indskel_maybeinside= C(1,indskel_maybeinside)'*z_axis(1)+C(2,indskel_maybeinside)'*z_axis(2)+C(3,indskel_maybeinside)'*z_axis(3);

% checking if the neighbouring points are above or below spherical cap
if strcmp(direction,'above')
    relind = zcoord_indskel_maybeinside>zcoord_proxy;
else
    relind = zcoord_indskel_maybeinside<zcoord_proxy;
end

indskel_inside = indskel_maybeinside(relind); % near field

% PLOTTING Gamma_near and Gamma_self:
plot3(C(1,indskel), C(2,indskel), C(3,indskel),'r.','MarkerSize',9);
plot3(C(1,indskel_inside), C(2,indskel_inside), C(3,indskel_inside), 'k.','MarkerSize',9);

text(0.7, -0.4, 1.1, '\Gamma_{near}', 'FontSize',17);
text(0.2, -0.25,1.2, '\Gamma_{self}','FontSize', 17,'Color','r');
text(0.56, 0.7, 0.6, '\Gamma_{far}', 'FontSize', 17,'Color','k');

% *** for exporting 
% addpath('export_fig-master'); 
% export_fig 'SourcesSelfNearFar' -png -m4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NODES = LOCAL_get_tree(C,nbox_max)

% NODES holds the binary tree for the recursive algorithm in a cell
% structure
% Blocks of the binary tree are ordered as in Section 3.4 (3.12)
% NODES is encoded as follows:

%--- TREE INFORMATION
%   NODES(01,ibox) = [pointer to geometry information]
%   NODES(02,ibox) = ilevel
%   NODES(03,ibox) = ifather
%   NODES(04,ibox) = list of children
%   NODES(05,ibox) = number of children
%   NODES(06,ibox) = index of first point in block
%   NODES(07,ibox) = nbox    - number of points in block
%   NODES(09,ibox) = kskel   - rank of off-diagonal block

%   NODES(11,ibox) = indices of skeleton points on contour, C
%   NODES(12,ibox) = matrix T output from the ID (holds linear combinations)
%   NODES(13,ibox) = indices of skeleton points in off-diagonal block
%

N = size(C,2);

NODES_0 = cell(13,N);

% Construct the top node.
NODES_0{1,1} = [1,N];
NODES_0{2,1} = 0;
NODES_0{3,1} = -1;
NODES_0{4,1} = [];
NODES_0{5,1} = 0;
NODES_0{6,1} = 1;
NODES_0{7,1} = N;
NODES_0{8,1} = -1;
NODES_0{9,1} = -1;

% Create smaller nodes via hierarchical subdivision.
% We sweep one level at a time.
ibox_last =  0;
ncreated  =  1;
while (ncreated > 0)
    ibox_first = ibox_last + 1;
    ibox_last  = ibox_last + ncreated;
    ncreated   = 0;
    for ibox = ibox_first:ibox_last
        nbox = NODES_0{7,ibox};
        if (nbox > nbox_max)
            nhalf              = floor(NODES_0{7,ibox}/2);
            if (nhalf > 0)
                ison1            = ibox_last + ncreated + 1;
                NODES_0{4,ibox}  = [NODES_0{4,ibox},ison1];   % Add a child to ibox.
                NODES_0{5,ibox}  = NODES_0{5,ibox}+1;         % Add a child to ibox.
                NODES_0{2,ison1} = NODES_0{2,ibox}+1;         % Set the level           of the child.
                NODES_0{3,ison1} = ibox;                      % Set the parent          of the child.
                NODES_0{4,ison1} = [];                        % Set the child list      of the child.
                NODES_0{5,ison1} = 0;                         % Set the child list      of the child.
                NODES_0{6,ison1} = NODES_0{6,ibox};           % Set the left-most index of the child.
                NODES_0{7,ison1} = nhalf;                     % Set the number of nodes of the child.
                ncreated         = ncreated + 1;
            end
            if (nhalf < NODES_0{7,ibox})
                ison2            = ibox_last + ncreated + 1;
                NODES_0{4,ibox}  = [NODES_0{4,ibox},ison2];   % Add a child to ibox.
                NODES_0{5,ibox}  = NODES_0{5,ibox}+1;         % Add a child to ibox.
                NODES_0{2,ison2} = NODES_0{2,ibox}+1;         % Set the level           of the child.
                NODES_0{3,ison2} = ibox;                      % Set the parent          of the child.
                NODES_0{4,ison2} = [];                        % Set the child list      of the child.
                NODES_0{5,ison2} = 0;                         % Set the child list      of the child.
                NODES_0{6,ison2} = NODES_0{6,ibox} + nhalf;   % Set the left-most index of the child.
                NODES_0{7,ison2} = NODES_0{7,ibox} - nhalf;   % Set the number of nodes of the child.
                ncreated = ncreated + 1;
            end
        end
    end
end

% This whole section is devoted to the simple task of
% decreasing the size of NODES_0.
% The operation NODES_0 = NODES_0{:,1:nboxes} apparently
% does not work for "cell" structures.
nboxes = ibox_last;
NODES  = cell(13,nboxes);
for ibox = 1:nboxes
    for j = 1:size(NODES_0,1)
        NODES{j,ibox} = NODES_0{j,ibox};
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xxc, max_dist] = LOCAL_get_circum_circle(Gamma_tau)
%--------------------------------------------------------------------------
% This function obtains the center and radius of the spherical cap which
% surrounds Gamma_tau
% Inputs:  Gamma_tau : the section of contour we would like to construct
%                      the circle around
% Outputs: xxc : the center of the circle
%          R   : the radius of the circle
%--------------------------------------------------------------------------

% Computing midpoint in 3D based on endpoints of Gamma_tau:
xxc = 0.5*(Gamma_tau([1,2,3],1) + Gamma_tau([1,2,3],end));
[th0,phi0,~]=cart2sph(xxc(1), xxc(2), xxc(3));

% Projecting onto sphere
[xxc(1), xxc(2),xxc(3)] = sph2cart(th0,phi0,1.0);
distsquare = (Gamma_tau(1,:)-xxc(1)).^2 + (Gamma_tau(2,:)-xxc(2)).^2+(Gamma_tau(3,:)-xxc(3)).^2;
max_dist = sqrt(max(distsquare));

% If the result is an absurdly large circle, then we instead base xxc on
% the center of mass of the points.
if ( (1.2*max_dist*max_dist) > ( (Gamma_tau(1,1) - Gamma_tau(1,end))^2 + (Gamma_tau(2,1) - Gamma_tau(2,end))^2 + (Gamma_tau(3,1)-Gamma_tau(3,end))^2))
    nloc = size(Gamma_tau,2);
    xxc = (1/nloc)*[sum(Gamma_tau(1,:)); sum(Gamma_tau(2,:)); sum(Gamma_tau(3,:))];
    [th0,phi0,~]=cart2sph(xxc(1), xxc(2), xxc(3));
    
    % Projecting onto sphere
    [xxc(1), xxc(2),xxc(3)] = sph2cart(th0,phi0,1.0);
    distsquare = (Gamma_tau(1,:)-xxc(1)).^2 + (Gamma_tau(2,:)-xxc(2)).^2+(Gamma_tau(3,:)-xxc(3)).^2;
    max_dist = sqrt(max(distsquare));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Cproxy] = LOCAL_construct_cap(xxc, dist, nproxy)
%--------------------------------------------------------------------------
% This function constructs the proxy spherical cap around Gammma_tau
% Inputs: xxc  : center of spherical cap
%         dist : distance between proxy cap and xxc
%       nproxy : number of points on proxy circle
%
tt = linspace(0, 2*pi*(1-1/nproxy), nproxy);

[th0,phi0,~]=cart2sph(xxc(1), xxc(2), xxc(3));
theta_c=acos(1-(dist^2)/2);
rproxy=sin(theta_c);

%r(t) % parametric equation for cap
r1=rproxy*cos(tt);
r2=rproxy*sin(tt);
r3=sqrt(1-rproxy^2);

% n(t) % outward normal
n1=sqrt(1-rproxy^2)*cos(tt);
n2=sqrt(1-rproxy^2)*sin(tt);
n3=-rproxy;

% s(t) % tangent vector
s1=-sin(tt);
s2=cos(tt);
s3=0;

% np(t) % principal normal
np1=-cos(tt);
np2=-sin(tt);
np3=zeros(1,length(tt));

% ds'(t) % arclength
ds=rproxy*ones(1,length(tt));

% curv   % curvature
curv=(1/rproxy)*ones(1,length(tt));

% forming new axes centred at xxc
[z1,z2,z3] = sph2cart(th0,phi0,1.0);
z_axis = [z1,z2,z3];
[x1,x2,x3] = sph2cart(th0,phi0-pi/2,1.0);
x_axis = [x1,x2,x3];
y_axis = cross(z_axis,x_axis);

% shifting cap to new origin on contour
Cproxy(1,:) = r1*x_axis(1) + r2*y_axis(1) + r3*z_axis(1);
Cproxy(2,:) = r1*x_axis(2) + r2*y_axis(2) + r3*z_axis(2);
Cproxy(3,:) = r1*x_axis(3) + r2*y_axis(3) + r3*z_axis(3);

Cproxy(7,:) = n1*x_axis(1) + n2*y_axis(1) + n3*z_axis(1);
Cproxy(8,:) = n1*x_axis(2) + n2*y_axis(2) + n3*z_axis(2);
Cproxy(9,:) = n1*x_axis(3) + n2*y_axis(3) + n3*z_axis(3);

Cproxy(10,:) = s1*x_axis(1) + s2*y_axis(1) + s3*z_axis(1);
Cproxy(11,:) = s1*x_axis(2) + s2*y_axis(2) + s3*z_axis(2);
Cproxy(12,:) = s1*x_axis(3) + s2*y_axis(3) + s3*z_axis(3);

Cproxy(16,:) = np1*x_axis(1) + np2*y_axis(1) + np3*z_axis(1);
Cproxy(17,:) = np1*x_axis(2) + np2*y_axis(2) + np3*z_axis(2);
Cproxy(18,:) = np1*x_axis(3) + np2*y_axis(3) + np3*z_axis(3);

Cproxy(19,:)=ds;
Cproxy(20,:)=curv;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LIST_NEI = LOCAL_construct_potential_neighborlist(C,NODES,dist_rel)

nboxes = size(NODES,2);

% Note that in the construction, we temporarily include
% a box in its own list of neighbors.
LIST_NEI    = cell(1,nboxes);
LIST_NEI{1} = 1;
for ibox = 2:nboxes
    [xxc, max_dist] = LOCAL_get_circum_circle(C(:,NODES{6,ibox}-1+(1:NODES{7,ibox})));
    % Collect of list of all the father's neighbors' sons:
    potneis = [];
    for jbox = LIST_NEI{NODES{3,ibox}};
        potneis = [potneis,NODES{4,jbox}];
    end
    % Loop over all potential neighbors to check which ones are actually "close".
    for jbox = potneis
        ind    = NODES{6,jbox} - 1 + (1:NODES{7,jbox});
        distsq = (C(1,ind)-xxc(1)).^2 + (C(2,ind)-xxc(2)).^2 + (C(3,ind)-xxc(3)).^2;
        if (min(distsq) < dist_rel*dist_rel*max_dist*max_dist)
            LIST_NEI{ibox} = [LIST_NEI{ibox},jbox];
        end
    end
end
% Finally we remove the box itself from the list of neighbors.
for ibox = 2:nboxes
    j = find(LIST_NEI{ibox} == ibox);
    LIST_NEI{ibox} = sort(LIST_NEI{ibox}([1:(j-1),(j+1):end]));
end

return