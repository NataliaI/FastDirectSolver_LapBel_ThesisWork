function [t_fmm, sigma, psi_hat] = FMMSolve(R, zeta, nbk, D, F, E, diagK,...
                                            dth, dzeta, Ck, Normal, dsda, acc,...
                                            x0, x1, gamma)
%--------------------------------------------------------------------------
% LAPLACE BELTRAMI SOLVER
% -Adapted from lapBel_driver.m and Laplace_Beltrami.m 
% -Uses GMRES to solve the Laplace-Beltrami equation on a simply connected 
% domain on the unit sphere, via an integral equation
% -Used for simply-connected domain, but can also be used for
% multiply-connected domain
 
% INPUTS: 
% R_FMM          : parametrized boundary curve (3 by N matrix holding
%                  coords) (output from island_geometry.m)
% zeta           : stereographic projection of boundary curve R
% nbk            : total number of points along all boundary contours
%                : in this case for a simply-connected contour nbk is the
%                : number of points along one contour, N
% D, F, E, diagK : blocks of system matrix obtained from functions 
%                  build_system.m and island_geometry.m
% dth            : grid size (h=2*pi/N) (obtained from island_geometry.m)
% dsda           : arc length 
% dzeta          : grid info for stereographic projection (output from
%                  island_geometry.m)
% Ck             : center of island in Cartesian coords
% Normal         : outward normal vector to boundary
% acc            : level of accuracy for GMRES
% x0             : point vortex location in Cartesian coords
% x1             : point vortex location on island in Cartesian coords
% gamma          : vortex strength

% OUTPUTS: 
% t_fmm          : time taken to solve system with FMM
% sigma          : solution to density function sigma
% psi_hat        : solution to the regular part of the stream function at
%                  x0
%--------------------------------------------------------------------------

nbod = 1; % number of islands 
f    = exact_solution_vort(R(1,:), R(2,:),R(3,:), x0, x1, gamma); % rhs for system
f(nbk+1: nbk+nbod) = zeros(nbod, 1);
 
% Solve via GMRES/FMM
Sinv = inv(D - F*E);  % Inverse of Schur complement for preconditioning
tic;
tol_gmres = acc;
soln_fmm = gmres(@(x) matvec(x, nbk, dth, zeta, dzeta, ...
                             diagK, E, F, D), f, nbk, tol_gmres, ...
                             nbk, @(x) leftPrec(x, nbk, E, F, Sinv));
t_fmm = toc;
sigma = soln_fmm(1: nbk)';   
  A_k = soln_fmm(nbk+1:end)';

% Evaluating solution 
u_calc  = double_layer_eval(dth, Ck, R, Normal, dsda, sigma, A_k, x0');
%figure; plot(sigma);
% psi_hat 
psi_hat = (1/2)*(u_calc...
          +(gamma/(2*pi))*log(sqrt((x0(1)-x1(1)).^2+(x0(2)-x1(2)).^2+(x0(3)-x1(3)).^2))...
          +(gamma/(4*pi))*log(2));
    
end