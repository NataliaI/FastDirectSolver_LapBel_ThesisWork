
function f = exact_solution_vort(x, y, z, x0, x1, gamma)
%--------------------------------------------------------------------------
% Boundary data for point vortex motion 
% for function u in psi_hat = gamma*G(x,x0) - gamma* G(x,x1) +u 
% INPUTS: 
% x,y,z - coordinates where right hand side is evaluated along
% x0    - point vortex location in Omega
% x1    - point vortex location on island
% gamma - vortex strength 
% OUTPUTS:
% f     - rhs for linear system [I + K] sigma = f
%--------------------------------------------------------------------------
f =  (gamma/(2*pi))*log(sqrt((x-x0(1)).^2+(y-x0(2)).^2+(z-x0(3)).^2))...
    -(gamma/(2*pi))*log(sqrt((x-x1(1)).^2+(y-x1(2)).^2+(z-x1(3)).^2));
f = f';
end

