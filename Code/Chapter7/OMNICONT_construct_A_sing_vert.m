function A = OMNICONT_construct_A_sing_vert(Cproxy,ind1,C2,ind2,h,N,centers)
%--------------------------------------------------------------------------
% Constructs reduced vertical off-diagonal block of A=I+(1/pi)K with proxy
% points
% Proxy points correspond to row indices, or targets.
%  Inputs :
%  Cproxy : curve info of proxy circle
%    ind1 : row indices of proxy circle (1:nproxy)
%      C2 : curve info of boundary
%    ind2 : column indices of boundary
%       h : grid size
%       N : total number of points on each island
% centers : centers of islands (where sources are placed)
% Outputs :
%       A : vertical off-diagonal matrix block
%--------------------------------------------------------------------------

% constructing kernel K(r_i, r_j)
%---------------------------------
M = size(centers,1);

i_tot = length(ind1);
j_tot = length(ind2);

[Y1, X1] = meshgrid(C2(1,ind2), Cproxy(1,ind1));
[Y2, X2] = meshgrid(C2(2,ind2), Cproxy(2,ind1));
[Y3, X3] = meshgrid(C2(3,ind2), Cproxy(3,ind1));

denom = (X1-Y1).^2+(X2-Y2).^2+(X3-Y3).^2;
n1 = repmat((C2(7,ind2)),i_tot,1);
n2 = repmat((C2(8,ind2)),i_tot,1);
n3 = repmat((C2(9,ind2)),i_tot,1);
dsY = repmat((C2(19,ind2)),i_tot,1);

A = ((h/pi)*((((X1-Y1).*n1+(X2-Y2).*n2+(X3-Y3).*n3))./(denom)).*dsY);

% adding on extra log sources
%-----------------------------
Sources = zeros(size(A));
for k = 2:M % going through islands and adding log sources
    indj   = 1:j_tot;
    % picking out column indices coresponding to island k
    indj_k = indj(ind2 >= (k-1)*N + 1 & ind2 <= k*N);
    ind2_k = ind2(ind2 >= (k-1)*N + 1 & ind2 <= k*N);
    
    if ~isempty(indj_k)  % if there are column indices corresponding to the kth island:
        % then add on sources
        Sources(1:i_tot,indj_k) = [-2*h*(1/(2*pi))*log(sqrt((Cproxy(1,ind1)-centers(1,1)).^2+(Cproxy(2,ind1)-centers(1,2)).^2+(Cproxy(3,ind1)-centers(1,3)).^2)/4)]'*C2(19,ind2_k)...
                                  +[2*h*(1/(2*pi))*log(sqrt((Cproxy(1,ind1)-centers(k,1)).^2+(Cproxy(2,ind1)-centers(k,2)).^2+(Cproxy(3,ind1)-centers(k,3)).^2)/4)]'*C2(19,ind2_k);
    end
end

A=A+Sources;

return