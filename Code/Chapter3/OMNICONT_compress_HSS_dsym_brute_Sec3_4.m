function [Atilde, NODES] = OMNICONT_compress_HSS_dsym_brute_Sec3_4(C,nbox_max,acc,v)
%--------------------------------------------------------------------------
% Brute Force Compression
% Plots skeleton points on each level 
%
% Modified from:
% http://amath.colorado.edu/faculty/martinss/2014_CBMS/codes.html
% Inputs:
% C        : contains boundary curve information
% nbox_max : specifies the maximum number of points allowable on the finest
%            level of the binary tree
% acc      : level of precision for the ID
% v        : specifies the view for plots 
%
% Outputs:
% Atilde : Compressed matrix, embedded in a sparse matrix (Chapter 3.7)
% NODES  : cell structure which stores information from each recursive level
%          in binary tree
%        : described in function LOCAL_get_tree below
%--------------------------------------------------------------------------

% Computing tree structure.
N          = size(C,2);
NODES      = LOCAL_get_tree(C,nbox_max);
nboxes     = size(NODES,2);
num_levels = log2(nboxes+1)-1;
levels     = 2.^(1:num_levels);

% figure handles: 
h_diag = [];  h_offdiag = []; 

% For storing indices of skeleton points: 
CurrentSkeleton = []; 
PrevSkeleton = []; 

% Compressing all blocks, going from smaller to larger:

% Obtaining indices of : Gamma_tau = indskel
%                        Gamma_tau_complement= ind_offd
for ibox = nboxes:(-1):2
    ind_offd = [1:(NODES{6,ibox}-1),(NODES{6,ibox}+NODES{7,ibox}):N];
    if (NODES{5,ibox}==0) % ibox has no sons.
        indskel = NODES{6,ibox} - 1 + (1:NODES{7,ibox});
    elseif (NODES{5,ibox}==1) % ibox has precisely one son
        ison1   = NODES{4,ibox}(1);
        indskel = NODES{11,ison1};
    elseif (NODES{5,ibox}==2) % ibox has precisely two sons
        ison1   = NODES{4,ibox}(1);
        ison2   = NODES{4,ibox}(2);
        indskel = [NODES{11,ison1},NODES{11,ison2}];
    end
    
    % building horizontal and vertical off-diagonal block
    A21 = OMNICONT_construct_A_offd(C,ind_offd,indskel);
    A12 = OMNICONT_construct_A_offd(C,indskel,ind_offd);
    
    % Computing skeletons.
    [T,J] = id_decomp([A21;A12'],acc);
    k     = size(T,1);
    
    % Recording skeletons:
    NODES{11,ibox} = indskel(J(1:k));
    NODES{12,ibox} = T;
    NODES{13,ibox} = J;
    NODES{ 9,ibox} = k;
    
    CurrentSkeleton = [CurrentSkeleton indskel(J(1:k))]; 
    
    figure(1); hold on 
    if ~isempty(h_diag)
        delete(h_diag); 
    end 
    if ~isempty(h_offdiag)
        delete(h_offdiag); 
    end 
    
    % plotting skeleton points for Gamma_tau 
    % (red is used for horizontal skeletons) 
    % (can change to blue for (identical) vertical skeleton)
    h_diag = plot3(C(1,indskel(J(1:k))), C(2,indskel(J(1:k))), C(3, indskel(J(1:k))),'r.','MarkerSize',13);
    
    % plotting old skeleton from previous level for Gamma_tau complement
    if isempty(PrevSkeleton) % if we are on the 1st level then plot all far field indices
        h_offdiag = plot3(C(1, ind_offd), C(2,ind_offd), C(3, ind_offd), 'k.', 'MarkerSize', 13); 
    else % pick out all of the skeleton indices from the previous level that are in the far field 
        ind_prevskel = ismember(ind_offd, PrevSkeleton);
        indfar = ind_offd(ind_prevskel);
    
        h_offdiag = plot3(C(1,indfar), C(2,indfar), C(3,indfar), 'k.','MarkerSize',13);
        view(v); alpha(0.8); 
    end
    keyboard;
    disp('click continue or type return');
    
    % *** for exporting: 
    % str = ['BFRecLev' num2str(num_levels-NODES{2,ibox}+1) 'HorSkel' num2str(ibox)]; 
    % addpath('export_fig-master');
    % export_fig(str, '-png', '-m4');
    
    % plotting complete skeleton on each level
    if ismember(ibox, levels)
        PrevSkeleton = CurrentSkeleton; 
        CurrentSkeleton = []; 
        
        [th, ph] = meshgrid(linspace(0, pi,100),linspace(0,2*pi,200)); 
        x = cos(ph).*sin(th);
        y = sin(ph).*sin(th); 
        z = cos(th);
        
        figure; hold on
        surf(x,y,z); shading interp; alpha(0.8); 
        hold on 
        cmap = gray(100);
        colormap(cmap(50:100,:));
        axis equal; axis off; 
        set(gcf, 'Color', 'w');
        view(v);
        plot3(C(1,PrevSkeleton),C(2,PrevSkeleton),C(3,PrevSkeleton), 'r.','MarkerSize',13);
        
        % for exporting: 
        % str = ['BFRecLev' num2str(num_levels-NODES{2,ibox}+1) 'HorSkelAll'];
        % export_fig(str, '-png', '-m4');
       
    end
end

% Constructing sparse matrix Atilde
Atilde = OMNICONT_construct_Atilde(C,NODES);

return

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
