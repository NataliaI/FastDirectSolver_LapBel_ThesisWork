function A = OMNICONT_construct_A_diag_MC(C,ind,h,N,centers)
%--------------------------------------------------------------------------
% Constructs diagonal block of A=I+(1/pi)K

%  Inputs :
%      C : curve info of boundary
%    ind : row and column indices of boundary
%       h : grid size
%       N : total number of points on each island
% centers : centers of islands (where sources are placed)
% Outputs :
%       A : diagonal matrix block
%--------------------------------------------------------------------------

% constructing kernel K(r_i, r_i)
%---------------------------------

M = size(centers,1);

nloc = length(ind);
[Y1, X1] = meshgrid(C(1,ind), C(1,ind));
[Y2, X2] = meshgrid(C(2,ind), C(2,ind));
[Y3, X3] = meshgrid(C(3,ind), C(3,ind));
denom    = (X1-Y1).^2+(X2-Y2).^2+(X3-Y3).^2;

n1  = repmat((C(7,ind)),nloc,1);
n2  = repmat((C(8,ind)),nloc,1);
n3  = repmat((C(9,ind)),nloc,1);
dsY = repmat((C(19,ind)),nloc,1) + eye(nloc);

A = ((h/pi)*((((X1-Y1).*n1+(X2-Y2).*n2+(X3-Y3).*n3))./(denom)).*dsY);

for i=1:nloc
    A(i,i) = (-h/(2*pi))*C(19,ind(i))*[C(10,ind(i)), C(11,ind(i)), C(12, ind(i))]*(C(20,ind(i))*[C(17, ind(i))*C(3,ind(i))-C(18,ind(i))*C(2,ind(i));...
                                       -(C(16,ind(i))*C(3,ind(i))-C(18,ind(i))*C(1,ind(i)));...
                                        C(16,ind(i))*C(2,ind(i))-C(17,ind(i))*C(1,ind(i))])+1;
end

% adding on extra log sources
%-----------------------------

Sources = zeros(size(A));

for k = 2:M % going through islands and adding log sources
    indj = 1:nloc;
    % picking out column indices coresponding to island k
    indj_k = indj(ind >= (k-1)*N+1 & ind <= k*N); % if corresponding to second contour
    ind_k  = ind(ind >= (k-1)*N+1 & ind <= k*N);
    
    if ~isempty(indj_k) % if there are column indices corresponding to the kth island:
        Sources(1:nloc,indj_k) = [-2*h*(1/(2*pi))*log(sqrt((C(1,ind)-centers(1,1)).^2+(C(2,ind)-centers(1,2)).^2+(C(3,ind)-centers(1,3)).^2)/4)]'*C(19,ind_k)...
                                 +[2*h*(1/(2*pi))*log(sqrt((C(1,ind)-centers(k,1)).^2+(C(2,ind)-centers(k,2)).^2+(C(3,ind)-centers(k,3)).^2)/4)]'*C(19,ind_k);
    end
end

A = A + Sources;

return