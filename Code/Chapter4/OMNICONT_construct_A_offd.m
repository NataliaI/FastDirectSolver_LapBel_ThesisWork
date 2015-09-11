function A = OMNICONT_construct_A_offd(C,ind1,ind2)
%--------------------------------------------------------------------------
% Constructs off-diagonal block of A=I+(1/pi)K
% Inputs:
% C:   curve info
% ind1: row indices of points along curve
% ind2: column indices of points along curve
%--------------------------------------------------------------------------

N     = size(C,2);
h     = 2*pi/N;

i_tot=length(ind1);
j_tot=length(ind2);

[Y1, X1]=meshgrid(C(1,ind2), C(1,ind1));
[Y2, X2]=meshgrid(C(2,ind2), C(2,ind1));
[Y3, X3]=meshgrid(C(3,ind2), C(3,ind1));

denom=(X1-Y1).^2+(X2-Y2).^2+(X3-Y3).^2;
n1=repmat((C(7,ind2)),i_tot,1);
n2=repmat((C(8,ind2)),i_tot,1);
n3=repmat((C(9,ind2)),i_tot,1);
dsY=repmat((C(19,ind2)),i_tot,1);

A=((h/pi)*((((X1-Y1).*n1+(X2-Y2).*n2+(X3-Y3).*n3))./(denom)).*dsY);

return