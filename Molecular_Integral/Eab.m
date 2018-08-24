function Eab = Eab(zeta,pos,l)
% amplitude of Hermite Gaussian on x, y, z axes

% input argument : 
% zeta = [a;b]
% pos = [A;B]
% L = [la;lb]
% input of function E
    ex = Eij(zeta,pos,l,1);
    ey = Eij(zeta,pos,l,2);
    ez = Eij(zeta,pos,l,3);
    [ex,ey,ez] = meshgrid(ey,ex,ez);
    Eab = ex.*ey.*ez;
end