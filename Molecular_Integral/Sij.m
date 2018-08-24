function S = Sij(zeta,pos,l,axes)
% Overlap Integral of 2 Gaussian
% input argument :
% zeta and pos are zeta and position of 2 Gaussian, respectively.
% l is angular quantum number, and axes is one of x,y,z coordinate(1, 2, 3, respectively)

E = Eij(zeta,pos,l,axes);
S = E(1)*sqrt(pi/sum(zeta,1));
end