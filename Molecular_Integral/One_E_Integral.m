function Output = One_E_Integral(zeta, L, pos_wave, pos_nuc, Z)
% one electron integral : 
% overlap integral, kinetic integral, and Coulomb integral

% input argument :
% Zeta is zeta of wave function, 
% pos_wave is the position of atom who own these wavefunction
% Z is charge of nuclie, pos_nuc is the position of nuclie

% output argument :
% One_E_Integral return a vector = 
% [kinetic integral, Coulomb integral, overlap integral]

p = sum(zeta,1); 
Pr = (zeta.'*pos_wave)./sum(zeta,1);
r_PC = pos_nuc-Pr;
% Overlap Integral
% Sab = S([i;j],1)*S([k;l],2)*S([m;n],3);
Sab =   Sij(zeta,pos_wave,L,1)...
       *Sij(zeta,pos_wave,L,2)...
       *Sij(zeta,pos_wave,L,3);
% Coulomb Integral
% sum(l,1) = la + lb
Clomb = sum(permute(sum(sum(Eab(zeta,pos_wave,L).*R0(sum(L,1),p,r_PC),1),2),[3 1 2]),1)...
    *(-2*pi*Z/p);
% Kinetic Integral
Tab = Tij(zeta,pos_wave,L,1)*Sij(zeta,pos_wave,L,2)*Sij(zeta,pos_wave,L,3)...
    + Sij(zeta,pos_wave,L,1)*Tij(zeta,pos_wave,L,2)*Sij(zeta,pos_wave,L,3)...
    + Sij(zeta,pos_wave,L,1)*Sij(zeta,pos_wave,L,2)*Tij(zeta,pos_wave,L,3);
% Output
Output = [Tab,Clomb,Sab];

end