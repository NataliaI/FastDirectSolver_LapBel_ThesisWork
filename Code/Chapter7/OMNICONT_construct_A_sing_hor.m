function A = OMNICONT_construct_A_sing_hor(C1,ind1,Cproxy,ind2,h,N,centers,indfar)
%--------------------------------------------------------------------------
% Constructs reduced horizontal off-diagonal block of A=I+(1/pi)K with proxy
% points
% Proxy points correspond to column indices, or sources.
%  Inputs :
%      C1 : curve info of boundary
%    ind1 : row indices of boundary
%  Cproxy : curve info of proxy circle
%    ind2 : column indices of proxy circle (1:nproxy)
%       h : grid size
%       N : total number of points on each island
% centers : centers of islands (where sources are placed)
%  indout : column indices of far-field (for determining which sources to
%           add)
% Outputs :
%       A : horizontal off-diagonal matrix block
%--------------------------------------------------------------------------

% constructing kernel K(r_i, r_j)
%---------------------------------
M = size(centers,1);

i_tot = length(ind1);
j_tot = length(ind2);

[Y1, X1] = meshgrid(Cproxy(1,ind2), C1(1,ind1));
[Y2, X2] = meshgrid(Cproxy(2,ind2), C1(2,ind1));
[Y3, X3] = meshgrid(Cproxy(3,ind2), C1(3,ind1));

denom = (X1-Y1).^2+(X2-Y2).^2+(X3-Y3).^2;
n1  = repmat((Cproxy(7,ind2)),i_tot,1);
n2  = repmat((Cproxy(8,ind2)),i_tot,1);
n3  = repmat((Cproxy(9,ind2)),i_tot,1);
dsY = repmat((Cproxy(19,ind2)),i_tot,1);

A = ((h/pi)*((((X1-Y1).*n1+(X2-Y2).*n2+(X3-Y3).*n3))./(denom)).*dsY);

% adding on extra log sources
%-----------------------------
Sources = zeros(size(A));

for k = 2:M  % going through islands and adding log sources
    % check indices of far field 
    indfar_k = indfar(indfar >=(k-1)*N + 1 & indfar <= k*N); 
    if ~isempty(indfar_k) % if colum indices correspond to contour k 
        % then add on sources
        Sources(1:i_tot,1:j_tot) = [-2*h*(1/(2*pi))*log(sqrt((C1(1,ind1)-centers(1,1)).^2+(C1(2,ind1)-centers(1,2)).^2+(C1(3,ind1)-centers(1,3)).^2)/4)]'*Cproxy(19,ind2)...
                                   +[2*h*(1/(2*pi))*log(sqrt((C1(1,ind1)-centers(k,1)).^2+(C1(2,ind1)-centers(k,2)).^2+(C1(3,ind1)-centers(k,3)).^2)/4)]'*Cproxy(19,ind2);
    end
end

A = A + Sources;
return