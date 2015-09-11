function A = OMNICONT_construct_A_sing(C1,ind1,C2,ind2,N)
%--------------------------------------------------------------------------
% Constructs reduced off-diagonal block of A=I+(1/pi)K with proxy points
% Inputs:
% C1, C2 : curve info - can be either original contour C or proxy circle
% ind1 : row indices of points along curve
% ind2 : column indices of points along curve
%    N : total number of points on original contour
%--------------------------------------------------------------------------

h  = 2*pi/N;

i_tot=length(ind1);
j_tot=length(ind2);

[Y1, X1]=meshgrid(C2(1,ind2), C1(1,ind1));
[Y2, X2]=meshgrid(C2(2,ind2), C1(2,ind1));
[Y3, X3]=meshgrid(C2(3,ind2), C1(3,ind1));

denom=(X1-Y1).^2+(X2-Y2).^2+(X3-Y3).^2;
n1=repmat((C2(7,ind2)),i_tot,1);
n2=repmat((C2(8,ind2)),i_tot,1);
n3=repmat((C2(9,ind2)),i_tot,1);
dsY=repmat((C2(19,ind2)),i_tot,1);

A=((h/pi)*((((X1-Y1).*n1+(X2-Y2).*n2+(X3-Y3).*n3))./(denom)).*dsY);

return