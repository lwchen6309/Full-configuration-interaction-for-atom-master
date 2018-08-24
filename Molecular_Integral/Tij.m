function T = Tij(zeta,pos,l,axes)
% Kinetic Integral of 2 Gaussian
% input argument :
% zeta and pos are zeta and position of 2 Gaussian, respectively.
% l is angular quantum number, and axes is one of x,y,z coordinate(1, 2, 3, respectively)

j = l(2,axes);
plusind = l;plusind(2,axes) = plusind(2,axes)+2;
minusind = l;minusind(2,axes) = minusind(2,axes)-2;
b = zeta(2);
T = -2*b^2*Sij(zeta,pos,plusind,axes)+b*(2*j+1)*Sij(zeta,pos,l,axes)...
    -j*(j-1)/2*Sij(zeta,pos,minusind,axes);
end