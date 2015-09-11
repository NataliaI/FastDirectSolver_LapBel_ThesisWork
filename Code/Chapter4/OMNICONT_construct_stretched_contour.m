function [C,t] = OMNICONT_construct_stretched_contour(flag_geom,params,contour_c,N, stretch)
%--------------------------------------------------------------------------
% The contour C is encoded as follows:
%  C(1,i), C(2,i), C(3,i)    -> r1 (t(i)), r2 (t(i)), r3 (t(i)) % counterclockwise
%  C(4,i), C(5,i), C(6,i)    -> r1'(t(i)), r2'(t(i)), r3'(t(i))
%  C(7,i), C(8,i), C(9,i)    -> n1 (t(i)), n2 (t(i)), n3 (t(i)) % inward normal
%  C(10,i), C(11,i), C(12,i) -> s1 (t(i)), s2 (t(i)), s3 (t(i)) % tangent (counterclockwise)
%  C(13,i), C(14,i), C(15,i) -> s1'(t(i)), s2'(t(i)), s3'(t(i))
%  C(16,i), C(17,i), C(18,i) -> np1(t(i)), np2(t(i)), np3(t(i)) % principal
%                                                                 normal
%  C(19,i)                   -> ds(t(i))
%  C(20,i)                   -> curv(t(i))                      % curvature
%
% ***** curve is taken counterclockwise so tangent is counterclockwise and 
% the normal is the inward normal 
%
% Inputs:
% flag_geom: flags the type of contour: 'star', 'circle', 'ellipse'
%
% params: contains parameters for different curves:
%   star   : [a w c]     a - degree of folding/curvature,
%                        w - number of arms
%                        c - size of star
%   cap : [r0 = sin(theta_c) theta_c]  theta_c - angle of elevation
%   ellipse: [a0 b0]   a0,b0 - major/minor axes
%   extended island: [L, delta]  L    - length of gap in island in stereographic
%                                        plane
%                                delta - (1+ delta) is the radius of circle in complex plane 
%
% contour_c : location of contour
% N : number of points on contour
%
% stretch : value by which to enlarge contour by (to avoid plotting close
%           to the boundary)
%
% Outputs:
% C : matrix of curve info
% t : discretized parameter t
% -------------------------------------------------------------------------

if strcmp(flag_geom, 'star')
    a = params(1); w=params(2); c=params(3);
    
    C = zeros(20,N);
    Ctemp=zeros(20,N);
    
    t = linspace(0,2*pi*(1 - 1/N),N);
    R = c*(1+a*cos(w*t));
    
    % r(t)   *** counterclockwise
    Ctemp(1,:) = R.*cos(t);
    Ctemp(2,:) = R.*sin(t);
    Ctemp(3,:) = sqrt(1-(Ctemp(1,:)).^2-(Ctemp(2,:)).^2);
    
    % r'(t)
    Rp     = -c*a*w*sin(w*t);
    Rpp    = -c*a*w^2*cos(w*t);
    Ctemp(4,:) = Rp.*cos(t)-R.*sin(t);
    Ctemp(5,:) = Rp.*sin(t)+R.*cos(t);
    Ctemp(6,:) = -R.*Rp./sqrt(1-R.^2);
    
    %s(t)  **** counterclockwise
    norm_rp = sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    Ctemp(10,:) = Ctemp(4,:)./norm_rp;
    Ctemp(11,:) = Ctemp(5,:)./norm_rp;
    Ctemp(12,:) = Ctemp(6,:)./norm_rp;
    
    % n(t)  *** inward normal 
    %n=cross([C(10,:)' C(11,:)' C(12,:)'],[C(1,:)' C(2,:)' C(3,:)']);
    n      = [Ctemp(11,:)'.*Ctemp(3,:)'-Ctemp(12,:)'.*Ctemp(2,:)'...
        -(Ctemp(10,:)'.*Ctemp(3,:)'-Ctemp(12,:)'.*Ctemp(1,:)')...
        Ctemp(10,:)'.*Ctemp(2,:)'-Ctemp(11,:)'.*Ctemp(1,:)'];
    norm_n = sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2);
    Ctemp(7,:) = (n(:,1)./norm_n)';
    Ctemp(8,:) = (n(:,2)./norm_n)';
    Ctemp(9,:) = (n(:,3)./norm_n)';
    
    %s'(t)
    r1pp    = Rpp.*cos(t)-2*Rp.*sin(t)-R.*cos(t);
    r2pp    = Rpp.*sin(t)+2*Rp.*cos(t)-R.*sin(t);
    r3pp    = -Rp.*Rp./sqrt(1-R.^2)-R.*(Rpp./sqrt(1-R.^2)+(Rp).^2.*R./((1-R.^2).^(3/2)));
    dnormp  = -(Ctemp(4,:).*r1pp+Ctemp(5,:).*r2pp+Ctemp(6,:).*r3pp)./((Ctemp(4,:).^2+Ctemp(5,:).^2+Ctemp(6,:).^2)).^(3/2);
    
    Ctemp(13,:) = r1pp./norm_rp+Ctemp(4,:).*dnormp;
    Ctemp(14,:) = r2pp./norm_rp+Ctemp(5,:).*dnormp;
    Ctemp(15,:) = r3pp./norm_rp+Ctemp(6,:).*dnormp;
    
    %np(t)
    norm_sp = sqrt((Ctemp(13,:)).^2+(Ctemp(14,:)).^2+(Ctemp(15,:)).^2);
    Ctemp(16,:) = Ctemp(13,:)./norm_sp;
    Ctemp(17,:) = Ctemp(14,:)./norm_sp;
    Ctemp(18,:) = Ctemp(15,:)./norm_sp;
    
    %ds'
    Ctemp(19,:) = sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    
    %curvature
    Ctemp(20,:) = norm_sp./norm_rp;
    
    C=TranslateContour(C,Ctemp, contour_c); % shifting contour to new origin
                                            % likely a better way to do
                                            % this...
    % Enlarging contour slightly:                                         
    C(1,:) = C(1,:)+stretch*C(7,:);
    C(2,:) = C(2,:)+stretch*C(8,:);
    C(3,:) = C(3,:)+stretch*C(9,:);
                                            
elseif strcmp(flag_geom, 'cap')
    C     = zeros(20,N);
    Ctemp = zeros(20,N);
    t     = linspace(0,2*pi*(1 - 1/N),N);
    
    r0=params(1);
    theta_c=params(2);
    
    % r(t)
    Ctemp(1,:)=r0*cos(t);
    Ctemp(2,:)=r0.*sin(t);
    Ctemp(3,:)=cos(theta_c)*ones(size(t));
    
    % r'(t)
    Ctemp(4,:)=-r0*sin(t);
    Ctemp(5,:)=r0*cos(t);
    Ctemp(6,:)=zeros(size(t));
    
    %s(t)
    norm_rp=sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    Ctemp(10,:)=Ctemp(4,:)./norm_rp;
    Ctemp(11,:)=Ctemp(5,:)./norm_rp;
    Ctemp(12,:)=Ctemp(6,:)./norm_rp;
    
    % n(t)
    %n=cross([C(10,:)' C(11,:)' C(12,:)'],[C(1,:)' C(2,:)' C(3,:)']);
    n=[Ctemp(11,:)'.*Ctemp(3,:)'-Ctemp(12,:)'.*Ctemp(2,:)'...
        -(Ctemp(10,:)'.*Ctemp(3,:)'-Ctemp(12,:)'.*Ctemp(1,:)')...
        Ctemp(10,:)'.*Ctemp(2,:)'-Ctemp(11,:)'.*Ctemp(1,:)'];
    norm_n=sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2);
    Ctemp(7,:)=(n(:,1)./norm_n)';
    Ctemp(8,:)=(n(:,2)./norm_n)';
    Ctemp(9,:)=(n(:,3)./norm_n)';
    
    %s'(t)
    r1pp=-r0*cos(t);
    r2pp=-r0*sin(t);
    r3pp=zeros(size(t));
    dnormp=-(Ctemp(4,:).*r1pp+Ctemp(5,:).*r2pp+Ctemp(6,:).*r3pp)./((Ctemp(4,:).^2+Ctemp(5,:).^2+Ctemp(6,:).^2)).^(3/2);
    
    Ctemp(13,:)=r1pp./norm_rp+Ctemp(4,:).*dnormp;
    Ctemp(14,:)=r2pp./norm_rp+Ctemp(5,:).*dnormp;
    Ctemp(15,:)=r3pp./norm_rp+Ctemp(6,:).*dnormp;
    
    %np(t)
    norm_sp=sqrt((Ctemp(13,:)).^2+(Ctemp(14,:)).^2+(Ctemp(15,:)).^2);
    Ctemp(16,:)=Ctemp(13,:)./norm_sp;
    Ctemp(17,:)=Ctemp(14,:)./norm_sp;
    Ctemp(18,:)=Ctemp(15,:)./norm_sp;
    
    %ds'
    Ctemp(19,:)=sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    
    %curvature
    Ctemp(20,:)=norm_sp./norm_rp;
    
    C=TranslateContour(C,Ctemp, contour_c); % shifting contour to new origin
    
    % Enlarging contour slightly:                                         
    C(1,:) = C(1,:)+stretch*C(7,:);
    C(2,:) = C(2,:)+stretch*C(8,:);
    C(3,:) = C(3,:)+stretch*C(9,:);
    
elseif strcmp(flag_geom, 'ellipse')
    C     = zeros(20,N);
    Ctemp = zeros(20,N);
    t     = linspace(0,2*pi*(1 - 1/N),N);
    
    a=params(1); b=params(2);
    
    % r(t)
    Ctemp(1,:)=a*cos(t);
    Ctemp(2,:)=b*sin(t);
    Ctemp(3,:)=sqrt(1-(Ctemp(1,:)).^2-(Ctemp(2,:)).^2);
    
    % r'(t)
    Ctemp(4,:)=-a*sin(t);
    Ctemp(5,:)=b*cos(t);
    Ctemp(6,:)=((a^2-b^2)*sin(t).*cos(t))./(sqrt(1-a^2*(cos(t)).^2-b^2*(sin(t)).^2));
    
    % useful terms
    c1=sqrt(-a^2*b^2+b^2*(cos(t)).^2+a^2*(sin(t)).^2);
    c2=c1./(Ctemp(3,:));
    c3=(-a^2*b^2+b^2*(cos(t)).^2+a^2*(sin(t)).^2);
    c4=1-(Ctemp(1,:)).^2-(Ctemp(2,:)).^2;
    
    % n(t)
    Ctemp(7,:)=-b*cos(t)*(a^2-1)./c1;
    Ctemp(8,:)=-a*sin(t)*(b^2-1)./c1;
    Ctemp(9,:)=-a*b./c2;
    
    %s(t)
    norm_rp=sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    Ctemp(10,:)=Ctemp(4,:)./norm_rp;
    Ctemp(11,:)=Ctemp(5,:)./norm_rp;
    Ctemp(12,:)=Ctemp(6,:)./norm_rp;
    
    %s'(t)
    Ctemp(13,:)=-(a*cos(t).*(((cos(t)).^4)*(-2*a^2*b^2+a^4+b^4)+((cos(t)).^2)*(2*a^4*b^2-2*a^4-2*a^2*b^4+2*a^2*b^2)-a^4*b^2+a^4+2*a^2*b^4-2*a^2*b^2-b^4+b^2))./(c3.^(3/2).*c4.^(1/2));
    Ctemp(14,:)=-(b*sin(t).*(((cos(t)).^4)*(-2*a^2*b^2+a^4+b^4)+((cos(t)).^2)*(2*a^4*b^2-2*a^4-2*a^2*b^4+2*a^2*b^2)+a^2*b^4-2*a^2*b^2+a^2))./(c3.^(3/2).*c4.^(1/2));
    Ctemp(15,:)=-(a^2-b^2)*(((cos(t)).^4)*(a^2-b^2)+((cos(t)).^2)*(2*a^2*b^2-2*a^2)-a^2*b^2+a^2)./(c3.^(3/2));
    
    %np(t)
    norm_sp=sqrt((Ctemp(13,:)).^2+(Ctemp(14,:)).^2+(Ctemp(15,:)).^2);
    Ctemp(16,:)=Ctemp(13,:)./norm_sp;
    Ctemp(17,:)=Ctemp(14,:)./norm_sp;
    Ctemp(18,:)=Ctemp(15,:)./norm_sp;
    
    %ds'
    Ctemp(19,:)=sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    
    %curvature
    Ctemp(20,:)=norm_sp./norm_rp;
    C=TranslateContour(C,Ctemp, contour_c);
    
    % Enlarging contour slightly:                                         
    C(1,:) = C(1,:)+stretch*C(7,:);
    C(2,:) = C(2,:)+stretch*C(8,:);
    C(3,:) = C(3,:)+stretch*C(9,:);
    
elseif strcmp(flag_geom, 'extended island')   
    C     = zeros(20,N);
    Ctemp = zeros(20,N);
    t     = linspace(0,2*pi*(1 - 1/N),N);
    
    L    = params(1); 
    delta = params(2);
    
    % constants
    c1  = 4*L*(1 + delta)*(2 + 2*delta + delta^2);
    c2  = (4 + 8*delta + 4*delta^2);
    c3  = 4*delta^2 + 4*delta^3 + delta^4 + 8*L^2*delta + 4*L^2*delta^2 + 4*L^2; 
    c4  = -4*L*(1 + delta)*(2 + delta)*delta;
    c5  = 4 + 8*delta + 4*delta^2;
    c6  = 4*delta^2 + 4*delta^3 + delta^4 + 8*L^2*delta + 4*L^2*delta^2 + 4*L^2; 
    c7  = (-4*delta - 2*delta^2 - 2); 
    c8  = (-2 - 4*delta - 6*delta^2 - 4*delta^3 - delta^4 + 8*L^2*delta + 4*L^2*delta^2 + 4*L^2);
    c9  = 4*delta + 2*delta^2 + 2; 
    c10 = 2 + 4*L^2 + 8*L^2*delta + 4*L^2*delta^2 + 4*delta + 6*delta^2 + 4*delta^3 + delta^4;
    
    % r(t)
    Ctemp(1,:) = 2*c1*cos(t)./(c2*cos(2*t)+c2+2*c3);
    Ctemp(2,:) = 2*c4*sin(t)./(c5*cos(2*t)+c5+2*c6);
    Ctemp(3,:) = (c7*cos(2*t)+c8)./(c9*cos(2*t)+c10);
    
    % r'(t)
    Ctemp(4,:) = -2*c1*sin(t)./(c2*cos(2*t)+c2+2*c3)+4*c1*cos(t).*c2.*sin(2*t)./(c2*cos(2*t)+c2+2*c3).^2;
    Ctemp(5,:) = 2*c4*cos(t)./(c5*cos(2*t)+c5+2*c6)+4*c4*sin(t)*c5.*sin(2*t)./(c5*cos(2*t)+c5+2*c6).^2;
    Ctemp(6,:) = -2*c7*sin(2*t)./(c9*cos(2*t)+c10)+(2*(c7*cos(2*t)+c8)).*c9.*sin(2*t)./(c9*cos(2*t)+c10).^2;
    Ctemp(4,:) = -C(4,:);
    Ctemp(5,:) = -C(5,:);
    Ctemp(6,:) = -C(6,:);
    
    % s(t)
    norm_rp = sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    Ctemp(10,:) = Ctemp(4,:)./norm_rp;
    Ctemp(11,:) = Ctemp(5,:)./norm_rp;
    Ctemp(12,:) = Ctemp(6,:)./norm_rp;
    
    % n(t)
    n = [Ctemp(11,:)'.*Ctemp(3,:)'-Ctemp(12,:)'.*Ctemp(2,:)' -(Ctemp(10,:)'.*Ctemp(3,:)'-Ctemp(12,:)'.*Ctemp(1,:)') Ctemp(10,:)'.*Ctemp(2,:)'-Ctemp(11,:)'.*Ctemp(1,:)'];
    norm_n = sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2);
    Ctemp(7,:) = (n(:,1)./norm_n)';
    Ctemp(8,:) = (n(:,2)./norm_n)';
    Ctemp(9,:) = (n(:,3)./norm_n)';
    
    % s'(t)
    r1pp = -2*c1*cos(t)./(c2*cos(2*t)+c2+2*c3)...
           -8*c1*sin(t).*c2.*sin(2*t)./(c2*cos(2*t)+c2+2*c3).^2 ...
           +16*c1*cos(t).*c2^2.*sin(2*t).^2./(c2*cos(2*t)+c2+2*c3).^3 ...
           +8*c1*cos(t).*c2.*cos(2*t)./(c2*cos(2*t)+c2+2*c3).^2;
     
    r2pp = -2*c4*sin(t)./(c5*cos(2*t)+c5+2*c6)...
           +8*c4*cos(t).*c5.*sin(2*t)./(c5*cos(2*t)+c5+2*c6).^2 ...
           +16*c4*sin(t).*c5^2.*sin(2*t).^2./(c5*cos(2*t)+c5+2*c6).^3 ...
           +8*c4*sin(t).*c5.*cos(2*t)./(c5*cos(2*t)+c5+2*c6).^2;
 
    r3pp = -4*c7*cos(2*t)./(c9*cos(2*t)+c10)...
           -8*c7*sin(2*t).^2*c9./(c9*cos(2*t)+c10).^2 ...
           +(8*(c7*cos(2*t)+c8)).*c9^2.*sin(2*t).^2./(c9*cos(2*t)+c10).^3 ...
           +(4*(c7*cos(2*t)+c8)).*c9.*cos(2*t)./(c9*cos(2*t)+c10).^2;
       
    dnormp = -(Ctemp(4,:).*r1pp+Ctemp(5,:).*r2pp+Ctemp(6,:).*r3pp)./((Ctemp(4,:).^2+Ctemp(5,:).^2+Ctemp(6,:).^2)).^(3/2);
    
    Ctemp(13,:) = r1pp./norm_rp + Ctemp(4,:).*dnormp;
    Ctemp(14,:) = r2pp./norm_rp + Ctemp(5,:).*dnormp;
    Ctemp(15,:) = r3pp./norm_rp + Ctemp(6,:).*dnormp;
    
    % np(t)
    norm_sp = sqrt((Ctemp(13,:)).^2+(Ctemp(14,:)).^2+(Ctemp(15,:)).^2);
    Ctemp(16,:) = Ctemp(13,:)./norm_sp;
    Ctemp(17,:) = Ctemp(14,:)./norm_sp;
    Ctemp(18,:) = Ctemp(15,:)./norm_sp;
    
    % ds'
    Ctemp(19,:) = sqrt((Ctemp(4,:)).^2+(Ctemp(5,:)).^2+(Ctemp(6,:)).^2);
    
    % curvature
    C(20,:) = norm_sp./norm_rp;
    
    C = TranslateContour(C, Ctemp, contour_c);
    
    % Enlarging contour slightly:                                         
    C(1,:) = C(1,:)+stretch*C(7,:);
    C(2,:) = C(2,:)+stretch*C(8,:);
    C(3,:) = C(3,:)+stretch*C(9,:);
      
end    
return

function C=TranslateContour(C,Ctemp,xxc)
[th0,phi0,~]=cart2sph(xxc(1), xxc(2), xxc(3));
z_axis=xxc;
[x1, x2, x3]=sph2cart(th0, phi0-pi/2,1);
x_axis=[x1 x2 x3];
y_axis=cross(z_axis,x_axis);

% shifting to new origin
% r(t)
C(1,:)=Ctemp(1,:)*x_axis(1) + Ctemp(2,:)*y_axis(1) + Ctemp(3,:)*z_axis(1);
C(2,:)=Ctemp(1,:)*x_axis(2) + Ctemp(2,:)*y_axis(2) + Ctemp(3,:)*z_axis(2);
C(3,:)=Ctemp(1,:)*x_axis(3) + Ctemp(2,:)*y_axis(3) + Ctemp(3,:)*z_axis(3);

% r'(t)
C(4,:)=Ctemp(4,:)*x_axis(1) + Ctemp(5,:)*y_axis(1) + Ctemp(6,:)*z_axis(1);
C(5,:)=Ctemp(4,:)*x_axis(2) + Ctemp(5,:)*y_axis(2) + Ctemp(6,:)*z_axis(2);
C(6,:)=Ctemp(4,:)*x_axis(3) + Ctemp(5,:)*y_axis(3) + Ctemp(6,:)*z_axis(3);

% s(t)
C(10,:)=Ctemp(10,:)*x_axis(1) + Ctemp(11,:)*y_axis(1) + Ctemp(12,:)*z_axis(1);
C(11,:)=Ctemp(10,:)*x_axis(2) + Ctemp(11,:)*y_axis(2) + Ctemp(12,:)*z_axis(2);
C(12,:)=Ctemp(10,:)*x_axis(3) + Ctemp(11,:)*y_axis(3) + Ctemp(12,:)*z_axis(3);

% n(t)
C(7,:)=Ctemp(7,:)*x_axis(1) + Ctemp(8,:)*y_axis(1) + Ctemp(9,:)*z_axis(1);
C(8,:)=Ctemp(7,:)*x_axis(2) + Ctemp(8,:)*y_axis(2) + Ctemp(9,:)*z_axis(2);
C(9,:)=Ctemp(7,:)*x_axis(3) + Ctemp(8,:)*y_axis(3) + Ctemp(9,:)*z_axis(3);

%s'(t)
C(13,:)=Ctemp(13,:)*x_axis(1) + Ctemp(14,:)*y_axis(1) + Ctemp(15,:)*z_axis(1);
C(14,:)=Ctemp(13,:)*x_axis(2) + Ctemp(14,:)*y_axis(2) + Ctemp(15,:)*z_axis(2);
C(15,:)=Ctemp(13,:)*x_axis(3) + Ctemp(14,:)*y_axis(3) + Ctemp(15,:)*z_axis(3);

%np(t)
C(16,:)=Ctemp(16,:)*x_axis(1) + Ctemp(17,:)*y_axis(1) + Ctemp(18,:)*z_axis(1);
C(17,:)=Ctemp(16,:)*x_axis(2) + Ctemp(17,:)*y_axis(2) + Ctemp(18,:)*z_axis(2);
C(18,:)=Ctemp(16,:)*x_axis(3) + Ctemp(17,:)*y_axis(3) + Ctemp(18,:)*z_axis(3);

%ds'
C(19,:)=sqrt(Ctemp(4,:).^2+Ctemp(5,:).^2+Ctemp(6,:).^2);

%curvature
C(20,:)=Ctemp(20,:);

return



