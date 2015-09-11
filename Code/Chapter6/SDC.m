function [sol,t,sol_eul] = SDC(x0, x1, gamma, C, T, J, Mp, npanels, N, boundary, btilde, L, U, P, Q, R) 
%--------------------------------------------------------------------------
% This function implements the Spectral Deferred Correction method to solve
% for the trajectory of a vortex (Section 5.3.2) : 
% dx0(t)/dt = F(t,x0(t)); x0(0) = ICx0
% 
% Inputs: x0 - initial point vortex location 
%         x1 - location of point vortex on island 
%      gamma - point vortex strength
%         C  - boundary contour information
%         T  - final time T
%        Mp  - number of points per panel
%         J  - number of correction steps
%    npanels - number of panels
%         N  - total number of points along C
%   boundary - boundary condition 
%     btilde - vector of zeros to specify size of RHS 
%  L,U,P,Q,R - factors from UMFPACk - for solving system 
%
% Outputs: 
%        sol - solution x0(t) - single vortex trajectory
%          t - vector of time steps from t=[0,T]; 
%    sol_eul - initial Forward Euler approximation 
%--------------------------------------------------------------------------

h = 2*pi/N; 

ICx0=x0(1);       ICy0=x0(2);       ICz0=x0(3);       % initial condition 
ICx0_eul = x0(1); ICy0_eul = x0(2); ICz0_eul = x0(3); 

M = Mp*npanels;
pwidth = T/npanels;     % panel width
a = 0; Tp = a+pwidth;   % panel interval [a, Tp]
sol = zeros(M,3);       % for storing solutions for each panel
t = zeros(M,1); 
                       
[tp]   = chebpts(Mp);      % grid for one panel on [-1,1]

dt=zeros(length(tp)-1,1);  % step sizes
for i=1:Mp-1
    dt(i)=tp(i+1)-tp(i);
end

for p=1:npanels
    
    x      = zeros(Mp,3);
    x_eul  = zeros(Mp,3);
    
    % Forward Euler (working on t in [-1,1])
    x(1,:)     = [ICx0 ICy0 ICz0];              % storing solution (values of x0(t))
    x_eul(1,:) = [ICx0_eul ICy0_eul ICz0_eul];
    
    sigma_mat = zeros(N,Mp); % stores solutions for sigma for each time step
    for step = 1:Mp-1
        sigma = LOCAL_get_sigma(boundary, C, x0, x1, btilde, L, U, P, Q, R);
        sigma_mat(:,step) = sigma;   
        
        vel = LOCAL_eval_F(((Tp-a)/2)*tp(step)+(Tp+a)/2, x0, x1, gamma, sigma, C, h); % velocity on original interval [a,T]
      
        x(step+1,:)      = x(step,:)+dt(step)*((Tp-a)/2)*vel;   % stretches step h to interval [a,T]
        x_eul(step+1, :) = x_eul(step,:)+dt(step)*((Tp-a)/2)*vel;
        
        x0 = x(step+1,:);  % on original interval 
    end
    sigma = LOCAL_get_sigma(boundary, C, x0, x1, btilde, L, U, P, Q, R);  % plug into matrix
    sigma_mat(:,step+1) = sigma;
    
    S = cumsummat(Mp);  % spectral integration matrix
    
    for j = 1:J   % correction steps
        % compute residual function (for each t)
        vel = LOCAL_eval_F_mat((Tp/2)*(1+tp), x, x1, gamma, sigma_mat, C, h);
        residual = S*((Tp-a)/2)*vel-x+[ICx0*ones(Mp,1) ICy0*ones(Mp,1) ICz0*ones(Mp,1)];
        
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
    end
    
    % t = [t; (Tp-a)/2*tp+(Tp+a)/2];
    t((p-1)*Mp+1: p*Mp, :) = (Tp-a)/2*tp+(Tp+a)/2;
    sol((p-1)*Mp+1: p*Mp,:) = x; 
    sol_eul((p-1)*Mp+1: p*Mp,:) = x_eul; 
    
    % initial condition for next panel
    ICx0     = x(end,1);     ICy0     = x(end,2);     ICz0     = x(end,3); 
    ICx0_eul = x_eul(end,1); ICy0_eul = x_eul(end,2); ICz0_eul = x_eul(end,3);
    
    a=Tp;
    Tp=a+pwidth;
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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








