function A = OMNICONT_construct_A_diag(C,ind)
%--------------------------------------------------------------------------
% Constructs diagonal block of A=I+(1/pi)K
% Inputs: 
% C:   curve info
% ind: indices of points along curve 
%--------------------------------------------------------------------------

N     = size(C,2);
h     = 2*pi/N;

N=length(ind);
[Y1, X1]=meshgrid(C(1,ind), C(1,ind));
[Y2, X2]=meshgrid(C(2,ind), C(2,ind));
[Y3, X3]=meshgrid(C(3,ind), C(3,ind));
denom=(X1-Y1).^2+(X2-Y2).^2+(X3-Y3).^2;
n1=repmat((C(7,ind)),N,1);
n2=repmat((C(8,ind)),N,1);
n3=repmat((C(9,ind)),N,1);
dsY=repmat((C(19,ind)),N,1)+eye(N);

A=((h/pi)*((((X1-Y1).*n1+(X2-Y2).*n2+(X3-Y3).*n3))./(denom)).*dsY);
for i=1:N
    %A(i,i)=(-h/(2*pi))*C(19,ind(i))*[C(10,ind(i)), C(11,ind(i)), C(12, ind(i))]*(cross(C(20,ind(i))*[C(16,ind(i)) C(17,ind(i)) C(18,ind(i))],[C(1,ind(i)) C(2,ind(i)) C(3,ind(i))]))'+1;
    A(i,i)=(-h/(2*pi))*C(19,ind(i))*[C(10,ind(i)), C(11,ind(i)), C(12, ind(i))]*(C(20,ind(i))*[C(17, ind(i))*C(3,ind(i))-C(18,ind(i))*C(2,ind(i)); -(C(16,ind(i))*C(3,ind(i))-C(18,ind(i))*C(1,ind(i))); C(16,ind(i))*C(2,ind(i))-C(17,ind(i))*C(1,ind(i))])+1;
end

return